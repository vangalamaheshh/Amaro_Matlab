function [prj,m,D,V,Q,totvar]=pca(dat,K)
% [prj,m,D,V,Q]=pca(dat,K)
%    calculates the Principal Components (PCA) for a dataset.
%    dat - data, rows in dat are considered the data points
%    K - number of components
%    returns:
%    prj - the projection on the first K components
%    m - center of mass vector
%    D - has the eigenvalues on the diagonal
%    V - has eigenvectors as columns
%    Q - if not empty has the orthonormal basis used to project the
%    centered data before diagonalizing the covariance matrix
% 
% Example: 
%    p=pca(rand(100,10),3);
%    plot3(p(:,1),p(:,2),p(:,3),'x');
% 
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%


m=mean(dat);
mdat=dat-repmat(m,size(dat,1),1);
totvar = sum(mdat(:).^2/(size(mdat,1)-1));

% more coordinates than points
Q=[];
if (size(mdat,2)>=size(mdat,1))
  Q = orth(mdat'); 
  ndat = mdat*Q;
else
  ndat=mdat;
end


cv = cov(ndat);
if (nargin==1)
  [V,D]=eig(cv);
  [evs,evi]=sort(-diag(D));
  V=V(:,evi);
  D=D(evi,evi);
else 
  opt.disp=0;
  if size(cv,1)>=K
    [V,D] = eigs(cv,K,'LM',opt);
  else
    prj=-1; D=[]; V=[]; 
    return;
  end
end

V=V.*repmat(sign(median(V)),size(V,1),1);

prj = ndat*V;
if ~isempty(Q)
%  V= V'*Q';
   V= Q*V;
end




