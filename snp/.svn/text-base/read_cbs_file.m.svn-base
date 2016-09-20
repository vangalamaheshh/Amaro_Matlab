function C=read_cbs_file(fname,Craw,headerlines,use_cols,use_sample_names,read_zones,use_gp_names)

if ~exist('use_cols','var') || isempty(use_cols)
  use_cols=1:size(Craw.dat,2);
end

if ~exist('use_sample_names','var') || isempty(use_sample_names)
    fid = fopen(fname);
    str = fgetl(fid);
    for k = 1:headerlines
        str = fgetl(fid);
    end
    fclose(fid)
    str = textscan(str,'%s\t');
    chars = str{1}{1};
    isnum = regexp(chars,'[^\d]');
    if ~isempty(isnum)
      use_sample_names = 1;
    else
      use_sample_names=0;
    end
    
end

if ~exist('read_zones','var') || isempty(read_zones)
  read_zones=0;
end

C=[];
if iscell(fname)
  already_matched=0;
  Craw_left=Craw;
  for i=1:length(fname)
    if exist('headerlines','var')
      Cs{i}=read_cbs_file(fname{i},Craw_left,headerlines);
    else
      Cs{i}=read_cbs_file(fname{i},Craw_left);
    end
    Craw_left=reorder_D_cols(Craw_left,setdiff(1:size(Craw_left.dat,2),1:size(Cs{i}.dat,2)));
  end
  C=unite_Ds(Cs,'cols');
else
  fid=fopen(fname,'r');
  
  if exist('headerlines','var')
    if use_sample_names
      if read_zones
        cbs=textscan(fid,'%s%d%f%f%d%f%f%f','headerLines',headerlines,'TreatAsEmpty','NA');
      else
        cbs=textscan(fid,'%s%d%f%f%d%f','headerLines',headerlines,'TreatAsEmpty','NA');
      end
      snames=cellstr(unique_keepord(strvcat(cbs{1}),'rows'));
      nscans=length(snames);
      if isempty(Craw)
        Craw_matched.sdesc=snames;
        use_cols=1:length(snames);
      else
        if use_gp_names
            new_names=regexprep(Craw.sdesc,'-','.');   % replace - with .
            new_names=regexprep(new_names,'^([0-9])','X$1'); % if starts with a digit change to Xdigit
        else
            new_names=Craw.sdesc;
        end
        [Mt,m1,m2]=match_string_sets(snames,new_names);
        %keyboard
        Craw_matched=reorder_D_cols(Craw,m2);
      end
    else
      if read_zones
        cbs=textscan(fid,'%d%d%f%f%d%f%f%f','headerLines',headerlines);
      else
        cbs=textscan(fid,'%d%d%f%f%d%f%f%f','headerLines',headerlines);
      end
      nscans=double(max(cbs{1}));
      Craw_matched=reorder_D_cols(Craw,1:min(nscans,size(Craw.dat,2)));
    end
  else
    error('no header lines');
  end
  fclose(fid);
  C=Craw_matched;
  if isfield(C,'dat')
    C.cbs=NaN*ones(size(C.dat));
  end
  for i=1:nscans
    if use_sample_names
      scidx=strmatch(new_names{i},cbs{1},'exact');
    else
      scidx=find(cbs{1}==min(cbs{1})+use_cols(i)-1);
    end
    scidx=scidx(find(~isnan(cbs{3}(scidx))));
    rlend=cumsum(double(cbs{5}(scidx)));
    rlst=[1; rlend(1:(end-1))+1];
    cbs_rl=[ rlst rlend cbs{6}(scidx)];
    if read_zones
      zonechr_rl=[ rlst rlend cbs{7}(scidx)];
      zonegen_rl=[ rlst rlend cbs{8}(scidx)];
    end
    if ~isfield(C,'cbs')
      C.cbs=nan(size(C.dat,1),nscans);
    end
    nonnan=find(~isnan(C.dat(:,i)));
    if length(nonnan)~=size(C.dat,1)
      % has NaNs
      if max(cbs_rl(:,2),[],1)~=length(nonnan)
        warning([ 'nans dont match: ' C.sdesc{i} ' : ' num2str([ max(cbs_rl(:,2),[],1) length(nonnan)])]);
        nonnan=nonnan(1:max(cbs_rl(:,2),[],1));
      end
      tmp=derunlength(cbs_rl,nonnan);
    else
      tmp=derunlength(cbs_rl);
    end
    disp([i length(tmp) nnz(isnan(tmp))]);
    C.cbs(1:length(tmp),i)=tmp;
    C.cbs_rl{i}=cbs_rl;
    if read_zones
      C.zonechr_rl{i}=zonechr_rl;
      C.zonegen_rl{i}=zonegen_rl;
    end
  end
  if ~isfield(C,'marker')
    C.marker=cellstr(num2str((1:size(C.cbs,1))'));
  end
  C.cbs=single(C.cbs);
end


