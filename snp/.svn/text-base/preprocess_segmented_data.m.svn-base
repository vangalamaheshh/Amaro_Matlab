function CL=preprocess_segmented_data(C,segments_file,should_log,is_hmm,smooth_sz,peak_sz,min_delta,less_memory,use_cols,subtract_median,cnv_file)

if ~exist('segments_file','var') || isempty(segments_file)
  segments_file='output';
end

if ~exist('min_delta','var') || isempty(min_delta)
  min_delta=0;
end

if ~exist('less_memory','var') || isempty(less_memory)
  less_memory=0;
end

if ~exist('should_log','var') || isempty(should_log)
  should_log=1;
end

if ~exist('is_hmm','var') || isempty(is_hmm)
  is_hmm=0;
end

if ~exist('smooth_sz','var') || isempty(smooth_sz)
  smooth_sz=3;
end

if ~exist('peak_sz','var') || isempty(peak_sz)
  peak_sz=0; % we do not want to use recover_peaks
end

if ~exist('use_cols','var') 
  use_cols=[];
end

if ~exist('subtract_median','var') 
  subtract_median=0;
end

if ~is_hmm
  CL=read_cbs_file(segments_file,C,0,use_cols);
  if exist('cnv_file','var') && ~isempty(cnv_file)
    CL=remove_cnv(CL,cnv_file);
  end
  if should_log
    if abs(mean(CL.dat(1:100:end)))>0.5
      disp('taking the log of the data');
      disp('using 0.1 as a lower cutoff - REMOVE ME')
      CL.dat(CL.dat<0.1)=0.1;
      CL.dat=log2(CL.dat)-1;
    end
    if abs(mean(CL.cbs(1:100:end)))>0.5
      disp('taking the log of the segmented data');
      disp('using 0.1 as a lower cutoff - REMOVE ME')
      CL.cbs(CL.cbs<0.1)=0.1;
      CL.cbs=log2(CL.cbs)-1;
    end
  end
  if smooth_sz>0
    CL=smooth_cbs(CL,smooth_sz);
  end
  if min_delta>0
    CL=join_close(CL,min_delta,0.3);
  end
  
  if peak_sz>0
    CL=recover_peaks(CL,1,peak_sz);
  else
    CL.cbs_fixed=CL.cbs;
  end  
  
  CL.dat=CL.cbs_fixed;
  CL=rmfield_if_exists(CL,{'cbs_fixed','cbs'});
  if (0)
    CL.raw=C.dat; % save raw data (not used)
    % save intermediate steps (not used)
    CL.sm1=CL3.cbs;
    CL.sm2=CL3_sm.cbs;
    CL.sm2j=CL3_sm_j.cbs;  
    CL.sm3=CL3_sm_pk.cbs_fixed;
  end
else
  CL=C;
  CL.dat=CL.hmm;
  CL=rmfield_if_exists(CL,'hmm');
end

% CL.smooth=CL.dat;
CL.chrn=chromosome2num(CL.chr);

if subtract_median
  disp('subtracting median of nonX values');
  CL.medians=median(CL.dat(CL.chrn~=23,:),1);
  CL.dat=CL.dat-repmat(CL.medians,size(CL.dat,1),1);
% FIXME: update cbs_rl and raw

%  for i=1:length(CL3.cbs_rl)
%    CL3.cbs_rl{i}(:,3)=CL3.cbs_rl{i}(:,3)-CL3.medians(i);
%  end
end



