function C = extract_from_fwb_and_collapse(patfile,targlist,categdir,fwbfiles,P)
% extract_from_fwb_and_collapse(patfile,targlist,categdir,fwbfiles,P)
%
% Mike Lawrence 2010-10-14

import org.broadinstitute.cga.tools.seq.FixedWidthBinary;

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'include_gsc',false);
P = impose_default_value(P,'method','standard');

% sample list
C.file.patfile = patfile;
C.sample = load_struct(patfile);
demand_field(C.sample,'name');
C.ns = slength(C.sample);
if length(fwbfiles)~=C.ns, error('length(fwbfiles)~=C.ns'); end
C.sample.fwbfiles = fwbfiles; clear fwbfiles;

% load categories list file
C.file.categdir = categdir;
C.file.categlist = [categdir '/categs.txt'];
C.cat = load_struct(C.file.categlist);
C.cat = make_numeric(C.cat,'num');
if C.cat.num(1)~=1, error('First category should be 1!'); end
C.ncat = slength(C.cat);
if C.cat.num(end)~=C.ncat, error('Category numbers should not be skipped'); end

% open categories FWB file
C.file.categfwb = [categdir '/all.fwb'];
try
  fprintf('Reading file %s\n',C.file.categfwb);
  catfwb = FixedWidthBinary(C.file.categfwb);
  catfwb.setNullVal(0);
catch me
  fprintf('Error reading %s\n',C.file.categfwb);
  error(me.message);
end

% load target file
C.file.targlist = targlist;
C.targ = load_target_file(targlist);
demand_fields(C.targ,{'gene'});
C.nt = slength(C.targ);

% collapse targets to genes
fprintf('Collapsing targets to genes\n');
[C.gene.name tmp C.targ.gidx] = unique(C.targ.gene);
C.ng = slength(C.gene);

% check input files
fprintf('Checking input files... ');
fwb = FixedWidthBinary;
for s=1:C.ns, if ~mod(s,10), fprintf('%d/%d ',s,C.ns); end
  if ~iscell(C.sample.fwbfiles{s}), C.sample.fwbfiles{s} = {C.sample.fwbfiles{s}}; end
  C.sample.n_fwbfiles(s) = length(C.sample.fwbfiles{s});
  for i=1:C.sample.n_fwbfiles(s)
    try
      fwb.open(C.sample.fwbfiles{s}{i});
      fwb.close();
    catch me
      fprintf('Error opening %s\n',C.sample.fwbfiles{s}{i});
      error(me.message);
end,end,end
fprintf('\nAll input files OK\n');

if strcmp(P.method,'standard')

  % allocate storage for coverage
  C.cov.sc = zeros(C.ns,C.ncat);
  C.cov.gc = zeros(C.ng,C.ncat);
  C.cov.ts = zeros(C.nt,C.ns);
  C.cov.gs = zeros(C.ng,C.ns);
  if P.include_gsc, C.cov.gsc = zeros(C.ng,C.ns,C.ncat); end

  for s=1:C.ns
    fprintf('Processing sample %d/%d: %s\n',s,C.ns,C.sample.name{s});
    fh = cell(length(C.sample.fwbfiles{s}),1);
    for f=1:C.sample.n_fwbfiles(s)
      fprintf('\tReading file %s\n',C.sample.fwbfiles{s}{f});
      fh{f} = FixedWidthBinary(C.sample.fwbfiles{s}{f});
      fh{f}.setNullVal(0);
    end
    tt=tic; step = round(C.nt/100);
    for t=1:C.nt, if ~mod(t,step), fprintf('%d/%d ',t,C.nt); toc(tt), end
      disp(t);
      g = C.targ.gidx(t);
      cov = fh{1}.get(C.targ.chr(t),C.targ.start(t),C.targ.end(t));
      if C.sample.n_fwbfiles(s)>1
        for f=2:C.sample.n_fwbfiles(s)
          cov = cov + fh{f}.get(C.targ.chr(t),C.targ.start(t),C.targ.end(t));
      end,end
      cat = catfwb.get(C.targ.chr(t),C.targ.start(t),C.targ.end(t));
      tot = 0;
      for i=1:length(cov)
        if cov(i)>0
          c = cat(i);
          if c>=1 && c<=C.ncat
            C.cov.sc(s,c) = C.cov.sc(s,c) + 1;
            C.cov.gc(g,c) = C.cov.gc(g,c) + 1;
            if P.include_gsc, C.cov.gsc(g,s,c) = C.cov.gsc(g,s,c) + 1; end
            tot=tot+1;
      end,end,end
      C.cov.ts(t,s) = C.cov.ts(t,s) + tot;
      C.cov.gs(g,s) = C.cov.gs(g,s) + tot;
    end, fprintf('\n');
    for f=1:C.sample.n_fwbfiles(s), fh{f}.close(); end  
  end

elseif strcmp(P.method,'quick')
  % only tabulates C.cov.gc and C.terr.gc

  % allocate storage for coverage
  C.cov.gc = zeros(C.ncat,C.ng);

  % allocate storage for territory
  C.terr.gc = zeros(C.ncat,C.ng);

  for s=1:C.ns
    fprintf('Processing sample %d/%d: %s\n',s,C.ns,C.sample.name{s});
    fh = cell(length(C.sample.fwbfiles{s}),1);
    for f=1:C.sample.n_fwbfiles(s)
      fprintf('\tReading file %s\n',C.sample.fwbfiles{s}{f});
      fh{f} = FixedWidthBinary(C.sample.fwbfiles{s}{f});
      fh{f}.setNullVal(0);
    end
    tt=tic; step = round(C.nt/100);
    for t=1:C.nt, if ~mod(t,step), fprintf('%d/%d ',t,C.nt); toc(tt), end
      g = C.targ.gidx(t);
      cov = fh{1}.get(C.targ.chr(t),C.targ.start(t),C.targ.end(t));
      if C.sample.n_fwbfiles(s)>1
        for f=2:C.sample.n_fwbfiles(s)
          cov = cov + fh{f}.get(C.targ.chr(t),C.targ.start(t),C.targ.end(t));
      end,end
      cat = catfwb.get(C.targ.chr(t),C.targ.start(t),C.targ.end(t));
      if s==1
        % tabulate territory
        if length(cat)>1
           C.terr.gc(:,g) = C.terr.gc(:,g) + histc(cat,1:C.ncat);
        elseif length(cat)==1  % (special case--histc returns a row vector instead of column vector)
          if cat>=1 && cat<=C.ncat
            C.terr.gc(cat,g) = C.terr.gc(cat,g) + 1;
          end
        end
      end
      % tabulate coverage
      cat(cov==0) = [];
      if length(cat)>1
        C.cov.gc(:,g) = C.cov.gc(:,g) + histc(cat,1:C.ncat);
      elseif length(cat)==1  % (special case--histc returns a row vector instead of column vector)
        if cat>=1 && cat<=C.ncat
          C.cov.gc(cat,g) = C.cov.gc(cat,g) + 1;
        end
      end
    end, fprintf('\n');
    for f=1:C.sample.n_fwbfiles(s), fh{f}.close(); end
  end
  
  C.cov.gc = C.cov.gc';
  C.terr.gc = C.terr.gc';
    
else
  error('unknown P.method');
end

catfwb.close();



