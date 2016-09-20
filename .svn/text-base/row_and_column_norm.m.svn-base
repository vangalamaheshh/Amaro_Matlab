function m=row_and_column_norm(m,niter,unitvec)
% start with row

if nargin<2
  niter=10;
end

if nargin<3
  unitvec=0;
end

for i=1:niter
  if unitvec
    m1=dna_norm(m)./sqrt(size(m,2)-1);
    m1=dna_norm(m1')'./sqrt(size(m1,1)-1);
  else
    m1=dna_norm(m); %./sqrt(size(m,2)-1);
    m1=dna_norm(m1')'; % ./sqrt(size(m1,1)-1);
  end
  d=m1-m;
  m=m1;
  verbose(['mean=' num2str(mean(d(:))) '; std=' num2str(std(d(:)))]);
end
