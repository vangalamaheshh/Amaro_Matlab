function [M,P,hlarge,snps]=correct_batch_effect(M,min_sz,bonf_pv_thresh,absolute_pv)

if iscell(M)
  clear M;
  global M;
end  
  
batch_supid=strmatch('BATCH',M.supacc,'exact');
control_supid=strmatch('CTRL',M.supacc,'exact');
nbatch=max(M.supdat(batch_supid,:));
h=histc(M.supdat(batch_supid,:),1:nbatch);
hlarge=find(h>=min_sz);
P=zeros(size(M.dat,1),length(hlarge));
DC=zeros(size(M.dat,1),length(hlarge));

b=zeros(length(hlarge),size(M.dat,2));
for i=1:length(hlarge)
  %%% FIXME: add inner loop on tissue type (include normals in each tissue type)
  [p,s]=differential_analysis(M,find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(batch_supid,:)~=hlarge(i)),...
                              struct('method','ttest_minvar', ...
                                     'minvar',(0.4)^2),0);
%  if nnz(M.supdat(control_supid,M.supdat(batch_supid,:)==hlarge(i)))>0
%    DC(:,i)=mean(M.dat(:,find(M.supdat(control_supid,M.supdat(batch_supid,:)==hlarge(i)))),2)-...
%            mean(M.dat(:,find(M.supdat(control_supid,M.supdat(batch_supid,:)~=hlarge(i)))),2);
%  else
%    DC(:,i)=NaN;
%  end
            
  [p,s]=differential_analysis(M,find(M.supdat(batch_supid,:)==hlarge(i)),find(M.supdat(batch_supid,:)~=hlarge(i)),...
                              struct('method','ttest_minvar', ...
                                     'minvar',(0.4)^2),0);
    
  P(:,i)=p;
  b(i,find(M.supdat(batch_supid,:)==hlarge(i)))=1;
end

if exist('absolute_pv','var') && absolute_pv~=0
  disp(['Finding SNPs with P-value < ' num2str(absolute_pv) ' for at least one batch']);
  [minv,mini]=min(P,[],2);
  snps=find(minv<absolute_pv); 
%  keyboard
else
  disp(['Finding SNPs with Bonferroni corrected P-value < ' num2str(bonf_pv_thresh) ' for at least one batch']);
  P2=min(P*size(P,1),1);
  [minv,mini]=min(P2,[],2);
  snps=find(minv<bonf_pv_thresh);
end

disp(['correcting overall ' num2str(length(snps)) ' SNPs']);
mini=mini(snps);
disp(['number of SNPs corrected in each (large) batch :']);
disp([ hlarge' histc(mini,1:length(hlarge))]);
B=b(mini,:);
nB=sum(B,2);
x=M.dat(snps,:);
m1=sum(x.*B,2)./nB;
m2=sum(x.*(1-B),2)./(size(B,2)-nB);
delta=(m1-m2)./m2;
disp([ mean(abs(delta)) std(abs(delta))]);
hist(100*delta,50);
M.dat(snps,:)=M.dat(snps,:)-B.*repmat((m1-m2),1,size(B,2));


            
