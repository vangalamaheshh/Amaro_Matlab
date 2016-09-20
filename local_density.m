function [rho,d,r,nbr]=local_density(D,dist_params,k)

dist=dna_distance(D.dat,dist_params);
n=size(dist,1);
ds=sort(dist);
r=ds(k+1,:);
nbr=sparse(ds<=repmat(r,n,1));
x=log(ds(2:end,:));
y=repmat((log(1:(n-1)))',1,n);

xm=cumsum(x);
nm=(1:(n-1))';
ym=cumsum(ds);

for i=1:n
  i
  for j=2:(n-1)
    p=polyfit(x(1:j,i),y(1:j,i),1);
    rho(j,i)=p(2);
    d(j,i)=p(1);
  end
end

% generate a dendrogram
% connect bnr points if both are above t, scan t
