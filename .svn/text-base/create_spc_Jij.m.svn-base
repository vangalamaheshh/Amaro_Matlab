function [Jij,a]=create_spc_Jij(dat,K,a,is_dist,dist_type)

N=size(dat,1);
if is_dist
  
else
  if exist('dist_type','var')
      [S,D]=find_knn(dat,K,'and',1,1,0,dist_type);
  else
      [S,D]=find_knn(dat,K,'and',1,1,0);
  end
  D=triu(D,1);
end
[xi,xj,dij]=find(D);
if ~exist('a','var') || isempty(a)
  a=mean(dij);
end
kav=nnz(S)/N;
Jij=sortrows([xi xj 1/kav*(1/sqrt(2*pi)/a)*(exp(-(dij/a).^2/2))]); 
%Jij=sortrows([xi xj exp(-(dij/a).^2/2)]); 


