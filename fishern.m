function [V,X,SB,SW,score]=fishern(cl);
% vectors are rows
% cl is a cell array of matrices

N=length(cl);
D=size(cl{1},2);
n=zeros(N,1);
p=zeros(N,1);

for i=1:N
  n(i)=size(cl{i},1);
end
p=n./sum(n);

SW=zeros(D,D);
M=zeros(N,D);
for i=1:N
  SW=SW+p(i)*cov(cl{i},1);
  M(i,:)=mean(cl{i});
end

M0=sum(repmat(p,1,D).*M);
SB=zeros(D,D);
for i=1:N
  SB=SB+p(i)*(M(i,:)-M0)'*(M(i,:)-M0);
end
warning off
J=inv(SW)*SB;
warning on
%keyboard
score=trace(J);
[V,X] = eig(J);
[dum,di]=sort(-abs(diag(X)));
V=V(:,di);
X=X(di,di);






