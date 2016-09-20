function [trueres,randres,permpv]=gsea_phen_perm_by_parts(dat,phen,test,sets,nperm,nparts,lsfdir)

if nargin <6
  nparts=1;
end

if nparts>1
  n_gsea_tests=4;
  nsets=length(sets);
  n_gsea_tests=4;
  p=get_parts(1:nperm,nparts);
  l=lsf(lsfdir);
  h=zeros(nparts,1);
  for i=1:nparts
    [l,h(i)]=bsub(l,{'trueres','randres','permpv'},'gsea_phen_perm',{dat,phen,test,sets,length(p{i})});
  end
  [l,res]=wait(l); % wait for all
  trueres=zeros(nsets,n_gsea_tests);
  randres=zeros(nperm,nsets,n_gsea_tests);
  
  for i=1:nparts
    trueres=res{h(i)}.trueres;
    randres(p{i},:,:)=res{h(i)}.randres;
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
else
  [trueres,randres,permpv]=gsea_phen_perm(dat,phen,test,sets,nperm);
end





