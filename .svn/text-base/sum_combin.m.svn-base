function y=sum_combin(n,x)
% how many ways to break n into x bins

tmp=nchoosek(1:(n+x-1),x-1);
tmp=[ zeros(size(tmp,1),1) tmp repmat(n+x,size(tmp,1),1)];
y=diff(tmp,1,2);
y=y-1;
