function [trueres,randres,permpv]=gsea_phen_perm(dat,phen,test,sets,nperm)

nsets=length(sets);
n_gsea_tests=4;
res=zeros(nperm,nsets,n_gsea_tests);

use_samples=~isnan(phen);
dat=dat(:,use_samples);
phen=phen(use_samples);

[ks,kspv,rs,rspv,true_order]=gsea_step(dat,phen,test,sets);
trueres(:,1)=ks;
trueres(:,2)=kspv;
trueres(:,3)=rs;
trueres(:,4)=rspv;

randres=[];
permpv=[];
if nperm>0
  ns=length(phen);
  for i=1:nperm
    disp(i);
    r=randperm(ns);
    [ks,kspv,rs,rspv,order]=gsea_step(dat,phen(r),test,sets);
    randres(i,:,1)=ks;
    randres(i,:,2)=kspv;
    randres(i,:,3)=rs;
    randres(i,:,4)=rspv;
  end
  
  direction_vec=[1 -1 -1 -1];
  % ks - higher is better
  % kspv - lower is better
  % rs - lower is better
  % rspv - lower is better
  
  rr=randres.*repmat(reshape(direction_vec,1,1,4), ...
                     size(randres,1),size(randres,2));
  tr=repmat(reshape(trueres.*repmat(direction_vec,size(trueres,1),1),...
                    1,size(trueres,1),size(trueres,2)),[nperm 1 1]);
  
  permpv=squeeze(sum(rr>tr,1))/nperm;
end




