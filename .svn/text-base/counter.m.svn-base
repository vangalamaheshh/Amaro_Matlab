function c=counter(maxval,start_at_zero)

if ~exist('start_at_zero','var')
  start_at_zero=0;
end

c=(1:maxval(end))';
for i=(length(maxval)-1):-1:1
  sz_c=size(c,2);
  tmp=repmat(1:maxval(i),size(c,1),1);
  c=[ tmp(:) repmat(c,maxval(i),1) ];
end
