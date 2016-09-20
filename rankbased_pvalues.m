function [gp,fwer,fpr]=rankbased_pvalues(S,SR)

SSR=sort(SR);
[SS,idx]=sort(S);
[tmp,invidx]=sort(idx);
gp=zeros(size(S,1),1);
fwer=gp;
n=size(SR,2);
fpr=gp;
for i=1:size(S,1)
  gp(i)=length(find(SSR(i,:)>=SS(i)))/n;
end
gp=gp(invidx);
fwer=find_tail(S,SSR(1,:)',1)/n;
fpr=find_tail(S,SR(:),1)/n/size(S,1);



