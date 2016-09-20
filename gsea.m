function [ks,kspv,rs,rspv,ks_core,rs_core]=gsea(order,sets)

nsets=length(sets);
[tmp,revorder]=sort(order);
for i=1:nsets
  ro_sets{i}=revorder(sets{i});
end

ks=zeros(length(sets),1);
kspv=zeros(length(sets),1);
ks_core=zeros(length(sets),1);

n=length(order);
% Kolmogorov-Smirnov
for i=1:nsets
  [h,kspv(i),ks(i),ks_core(i)]=my_kstest2(ro_sets{i},1:n,[],1);
end

rs=zeros(length(sets),1);
rspv=zeros(length(sets),1);
rs_core=zeros(length(sets),1);

% Ranksum test
if nargout>2
  rs=zeros(length(sets),1);
  rspv=zeros(length(sets),1);
  for i=1:nsets
    [p,h,rsstruct]=my_ranksum(ro_sets{i},1:n);
    rs(i)=rsstruct.ranksum;
    rspv(i)=rsstruct.p1;
    rs_core(i)=-1;
  end
end






