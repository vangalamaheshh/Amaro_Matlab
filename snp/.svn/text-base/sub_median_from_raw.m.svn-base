function C=sub_median_from_raw(C,fix_cbs_rl)

if ~exist('fix_cbs_rl','var')
  fix_cbs_rl=1;
end

if mean(C.raw(:,1))>1
  disp('taking log of C.raw');
  C.raw(C.raw<0.1)=0.1; % FIXME: carry with C the lower cut-off (0.1)
  C.raw=log2(C.raw)-1;
end

% update C.medians or make they match C.cbs_rl
if isfield(C,'medians')
  med=C.medians;
  tmp=derunlength(C.cbs_rl);
  med2=median(tmp(C.chrn~=23,:));
  if any(med~=med2)
    disp('medians do not match');
  end
else
  disp('calculating medians');
  tmp=derunlength(C.cbs_rl);
  med=median(tmp);
  C.medians=med;
end

C.raw=C.raw-repmat(C.medians,size(C.raw,1),1);
if fix_cbs_rl
  for i=1:size(C.dat,2)
    C.cbs_rl{i}(:,3)=C.cbs_rl{i}(:,3)-C.medians(i);
  end
end
