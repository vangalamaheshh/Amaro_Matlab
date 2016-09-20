function [gi,gm,model]=find_ntop_genes(dat,n)

% data is log ratio to control

gm=mean(dat,2);
%cdat=dat-repmat(gm,1,size(dat,2));
cdat=dat;

gi=[];
for k=1:n
  
  dist_m=zeros(size(cdat,1),size(cdat,2));
  for r=1:size(cdat,2)
    r
    cdat_loo=cdat(:,setdiff(1:size(cdat,2),r));
    [prj,m,D,V,Q]=pca(cdat_loo,1);
    ncdat=dna_norm(cdat_loo);
    dist_m(:,r)=sum((ncdat-repmat(dna_norm(V'),size(ncdat,1), ...
                                  1)).^2,2);
  end
  [md,mi]=min(mean(dist_m,2));
  gi=[gi mi];
  disp(md)
  
  
end
