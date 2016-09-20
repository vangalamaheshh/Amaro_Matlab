function s=calc_copy_quality_score(C,quality_type,use_raw);

if ~exist('quality_type','var')
  quality_type='mean';
end

if ischar(quality_type)
  tmp.method=quality_type;
  quality_type=tmp;
end

if exist('use_raw','var') && use_raw
 % x=diff(C.raw,1,1);  
  C = calcdiff(C,'raw',1,1,'diff');  %put x in C as C.diff
else
   
    
   C = calcdiff(C,'dat',1,1,'diff');   %put x in C as C.diff
%  x=diff(C.dat,1,1);
  
  
end

switch quality_type.method
 case 'mean'
     C = itrfcn2(C,'diff','diffsquared',1,@(x) x.^2);
     s = itrfcn1(C,'diffsquared',1,@mean,1);
     C = rmfield(C,'diff');
     C = rmfield(C,'diffsquared');
%     s = mean(x.^2,1);
 case 'median'
     C = itrfcn2(C,'diff','diffsquared',1,@(x) x.^2);
     s = itrfcn1(C,'diffsquared',1,@median,1);
     C = rmfield(C,'diff');
     C = rmfield(C,'diffsquared');
%  s=median(x.^2,1);
 case 'medianabs'
     s = itrfcn1(C,'diff',1,@(x) nanmedian(abs(x),1)/(sqrt(2)*norminv(0.75)));
     C = rmfield(C,'diff');
%  s=nanmedian(abs(x),1)/(sqrt(2)*norminv(0.75));  
 case 'mad'
     
  x=[C.diff; -C.diff];
  s=mad(x,1,1)*mad_factor;
  C = rmfield(C,'diff');
 otherwise
  error('no such type');
end

