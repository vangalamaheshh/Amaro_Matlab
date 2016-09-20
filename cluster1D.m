function [G]=cluster1D(x,xc)
% x=rand(100,1); xc=0.01
n=0;
N=length(x);
[x1,m]=sort(x(:));
xL=x1-xc;
xH=x1+xc;
% low window edge
[q,iL]=histc(xL,x1); iL(iL<1)=1;
% high window edge
[q,iH]=histc(xH,x1); iH(iH<1)=length(x1);
% point density
xD=iH-iL;
class=1:N;
nxt=class;  % self-linked list
for i=1:N
    [xx,j] = max(xD(iL(i):iH(i))); j=j+iL(i)-1;
    [class,nxt]=connect1D(i,j,class,nxt);
end
GU=unique(class);
NG=length(GU);
G=0*class;
for i=1:NG
    G(class==GU(i))=i;
end
G(m)=G;


function [class,nxt]=connect1D(i,j,class,nxt)
% Saul & Knuth Equivalence Class algorithm
if (class(i)==class(j)), return;end
j1=j;
class(j1)=class(i);
while (nxt(j1)~=j)
    j1=nxt(j1);
    class(j1)=class(i);
end
i1=nxt(i);
nxt(i)=j;
nxt(j1)=i1;

