function [G]=clusterNN(xy,xyc)
%   ClusterNNN2D(xy,xyc)  
%   cluster in two dimensions of xy: x(:,1) 
%   with neighbor threshold xyc(1), and  xy(:,2) with threshold xyc(2)
%   inputs:
%       xy = xy(:,1) and xy(:,2) variables to be clustered
%       xyc = xc(1) and xc(2) cluster scale neighbor threshold 
%   outputs: 
%       G = group id for each cluster in order of x
%
% Author: Chip Stewart & Amaro Taylor-Weiner
% Created: 2013-08-30
%


G=[];
if isempty(xy), return; end
[N,dim]=size(xy);
if (dim<2), 
    [G]=clusterNN1d(xy,xyc)
    return; 
end
if (N<2), return; end
[Nc,dimc]=size(xyc);
if (Nc==1)
    xyc=repmat(xyc,N,1);
end
if (dimc==1)
    xyc=repmat(xyc,1,2);
end    
[Nc,dimc]=size(xyc);
if (Nc~=N)|(dimc~=dim)
    error('xyc dimensions not matched to xy')
    return
end


%cluster first dimension
%x1c=xc(1)
[G1]=clusterNN1d(xy(:,1),xyc(:,1));

% cluster second dimension
%x2c=xc(2)
[G2]=clusterNN1d(xy(:,2),xyc(:,2));

% combine groups
Gx = 10^(1+ceil(log10(max(G1))));
G0 = G2*Gx+G1;

% remap to consecutive group id's 
UG=unique(G0);
G=zeros(N,1);
for n=1:length(UG)
    k=find(G0==UG(n));
    G(k)=n+0*k;
end

if (0)
    UG=unique(G);
    for i=1:length(UG)
        fprintf(1,'%d ',i)
        k=find(G==UG(i));
        fprintf(1,'%.3f ',G(k))
        fprintf(1,'\n')
    end
end

function [G]=clusterNN1d(x,xc)
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
    [class,nxt]=connect(i,j,class,nxt);
end
GU=unique(class);
NG=length(GU);
G=0*class;
for i=1:NG
    G(class==GU(i))=i;
end
G(m)=G;

function [class,nxt]=connect(i,j,class,nxt)
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

function test
%%
clear
clf
%f='/Users/stewart/Projects/matlab/CCF_data_for_2dnn.mat'
%load(f)
%x=str2double(CCF_test.x)
%y=str2double(CCF_test.y)
%xc=(str2double(CCF_test.x_highccf)-str2double(CCF_test.x_lowccf))/3

for i=1:4

    
    subplot(2,2,i);
    N=round(10^(i/2));
    
    P0=[1 1; 1 0; 0.25 0.25; 0.5 0];
    C0=0.1;
    NG=size(P0,1);
    
    P=repmat(P0,N,1);
    x=randn(N*NG,1)*C0+P(:,1);
    xc=0*x+C0;
    y=randn(N*NG,1)*C0+P(:,2);
    yc=0*y+C0;
    
    plot(x,y,'o')
    G=clusterNN([x y],[xc yc]*1.1);
    ug=unique(sort(G));
    cm=hsv(length(ug))/2;
    for g=1:length(ug)
        k=find(G==ug(g));
        hold on;
        plot(x(k),y(k),'o','markeredgecolor','none','markerfacecolor',cm(g,:));
        text(0.1,0.9,sprintf('NG=%d, N=%d',length(ug),length(x)),'units','normalized');
    end
    plot(P(:,1),P(:,2),'kx','markersize',15,'linewidth',2);
    hold off
end

