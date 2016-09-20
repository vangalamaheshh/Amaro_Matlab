function [res,rep]=calc_reproducibility(D,nsplits_per_skewness,nperm,use_confounder,direction)
% assume first supdat is the phenotype and second is the confounder
addpath ~/projects/confounders

selection=struct('method','bhfdr','thresh',0.2);
if use_confounder
%  test_type=struct('method','ttest','nparts_perm',1,'nparts_fix', ...
%                   1,'confounding',D.supdat(2,:));
   test_type=struct('method','ftest', ...
                    'confounding',[]);
else
%  test_type=struct('method','ttest','nparts_perm',1,'nparts_fix', ...
%                   1);
   test_type=struct('method','ttest');
end

[splits,dw]=get_splits(D,nsplits_per_skewness,direction);
ns=length(splits);
res=cell(ns,2);
rep=zeros(ns,2);
rep(:,1)=1; %dw; % fixme!!!
for i=1:length(splits)
  % D1
  r.sp1=splits{i};
  D1=reorder_D_cols(D,r.sp1);
  if use_confounder
    test_type.confounding=D1.supdat(2,:);
  end  
  [r.idx,r.q,r.p,r.s,r.pi0,r.F]=get_top_markers(D1,1,test_type,nperm,selection);
  res{i,1}=r;
  % D2
  r.sp2=setdiff(1:size(D.dat,2),splits{i});
  D2=reorder_D_cols(D,r.sp2);
  if use_confounder
    test_type.confounding=D2.supdat(2,:);
  end  
  [r.idx,r.q,r.p,r.s,r.pi0,r.F]=get_top_markers(D2,1,test_type,nperm,selection);
  res{i,2}=r;
  
  rep(i,2)=length(intersect(res{i,1}.idx,res{i,2}.idx))/ ...
           length(union(res{i,1}.idx,res{i,2}.idx));
  v1=zeros(size(D.dat,1),1); v1(res{i,1}.idx)=1;
  v2=zeros(size(D.dat,1),1); v2(res{i,2}.idx)=1;
  [ct, chi2, pv]=crosstab(v1,v2);
  rep(i,3)=pv;
  disp([num2str(i) '/' num2str(length(splits)) ]);
  ct1=crosstab(D1.supdat(1,:),D1.supdat(2,:));
  ct2=crosstab(D2.supdat(1,:),D2.supdat(2,:));
  disp([ ct1 ct2 ])  
end


