function C=calc_segment_std_err(C,check_med)

if ~exist('check_med','var')
  check_med=0;
end

C.rl_std=cell(size(C.cbs_rl));
for i=1:size(C.dat,2)
  fprintf(1,'.');
  r=C.cbs_rl{i};
  s=zeros(size(r,1),1);
  for j=1:size(r,1)
    s(j)=mad(C.raw(r(j,1):r(j,2),i),1)*mad_factor/sqrt(r(j,2)-r(j,1)+1);
    if check_med
      tmp=median(C.raw(r(j,1):r(j,2),i),1);
      if tmp~=r(j,3)
        disp([i j tmp r(j,3)]);
%         error('segment does not match');
      end
    end
  end
  C.rl_std{i}=s;
end
disp('segment std err: Done.');
