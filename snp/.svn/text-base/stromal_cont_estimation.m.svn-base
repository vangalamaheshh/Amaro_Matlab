function [p,alpha,s2,cn,pval]=stromal_cont_estimation(r,ei,ce,delta)

if ~exist('delta','var')
  delta=[-1 1];
end

c=[ce-1 ce+1];
if (ei+delta(1))>=1
  beta(1)=r(ei+delta(1))/r(ei);
  pval(1)=(2*beta(1)-2)/((c(1)-2)-(ce-2)*beta(1));
else
  pval(1)=NaN;
  beta(1)=0;
end

if (ei+delta(2))<=length(r)
  beta(2)=r(ei+delta(2))/r(ei);
  pval(2)=(2*beta(2)-2)/((c(2)-2)-(ce-2)*beta(2));    
else
  pval(2)=NaN;
  beta(2)=0;
end
    
p=nanmean(pval,2);
alpha=(p*(ce-2)+2)/r(ei);
cn=(alpha.*r-2)/p+2;
s2=0;

if ei>=2
  s2=s2+((beta(1)*((ce-2)*p+2)-2)/p+2-c(1)).^2;
end
if ei<length(r)
  s2=s2+((beta(2)*((ce-2)*p+2)-2)/p+2-c(2)).^2;      
end





