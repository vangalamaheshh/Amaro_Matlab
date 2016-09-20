classdef bam
  properties
    chr
    st
    en
    R
    B
    S
    isz   % insert-size distribution
    fhm   % forward half-mapped reads
    rhm   % reverse half-mapped reads
    fhm_sm  % smoothed
    rhm_sm
    dhm   % difference between forward and reverse
    dhmz  % divided by constant
    shm   % sidedness
    shmz  % divided by constant
  end

  methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function obj = loadbam(obj,fname)
      [obj.R obj.B obj.S] = pull_from_bam(fname,obj.chr,obj.st,obj.en);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function obj = loadisz(obj,fname)
      tmp = load_matrix(fname);
      obj.isz = sum(tmp(:,2:end),2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate metrics:
% 1. number of forward reads starting here that have unmapped pairmate
% 2. number of backward reads ...

    function obj = calc(obj,sw,sw2,flag)
      if ~exist('sw','var'),sw = 300;end
      if ~exist('sw2','var'),sw2 = 50;end
      if ~exist('flag','var'), flag=true; end
      z = zeros(obj.en-obj.st+1,1);
      obj.fhm = z; obj.rhm = z; obj.dhm = z;
      for i=1:size(obj.R,1)
        % if read is paired and pairmate is unmapped
        if obj.R(i,3)>0 && obj.R(i,10)==-1
          pos = obj.R(i,4) - obj.st + 1;
          if pos>=1   % (if read starts in the window loaded)
            if obj.R(i,6)==0
              obj.fhm(pos) = obj.fhm(pos) + 1;
            elseif obj.R(i,6)==1
              obj.rhm(pos) = obj.rhm(pos) + 1;
            else error('illegal strand %d for mapped read\n',obj.R(i,6));
      end,end,end,end
      if flag
        obj.fhm_sm = smooth(obj.fhm,sw);
        obj.rhm_sm = smooth(obj.rhm,sw);
      end
      % forward: calculate at each point: the SUM in a sw ending at that point
      fcs = cumsum(obj.fhm);
      fws = zeros(length(fcs),1);
      fws(1+sw:end) = fcs(1+sw:end) - fcs(1:end-sw);
      % reverse: calculate at each point: the SUM in a sw beginning at that point
      rcs = cumsum(obj.rhm);
      rws = zeros(length(rcs),1);
      rws(1:end-sw) = rcs(1+sw:end) - rcs(1:end-sw);
      % combined: calculate at each point, dhm = [f(-) + r(+)] - [r(-) + f(+)]
      dhm = zeros(length(fcs),1);   % Note: having nan's slows down smooth ~100x
      dhm(1+sw:end-sw) = fws(1+sw:end-sw) + rws(1+sw:end-sw) - ...
                         fws(1+2*sw:end) - rws(1:end-2*sw);
      dhm = smooth(dhm,sw2);
      obj.dhm = dhm;
      % divide by median
      %med = median(abs(obj.dhm));
      %obj.dhmz = obj.dhm / med;
      % divide by constant
      obj.dhmz = obj.dhm / 5;
      % sidedness: calculate at each point, shm = [f(-) - f(+)] - [r(+) - r(-)]
      shm = zeros(length(fcs),1);   % Note: having nan's slows down smooth ~100x
      shm(1+sw:end-sw) = fws(1+sw:end-sw) - fws(1+2*sw:end) - ...
                         rws(1+sw:end-sw) + rws(1:end-2*sw);
      shm = smooth(shm,sw2);
      obj.shm = shm;
      obj.shmz = obj.shm / 5;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end
end
