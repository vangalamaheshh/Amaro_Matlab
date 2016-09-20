function [A,B,M,Z,s]=calc_delta(C,N,rng);
% go to range not just 1
% C is in raw scale
if ~exist('rng','var')
  rng=[0.05 1];
end
stp=(rng(2)-rng(1))/(N-1);

s=inf(size(C.dat,2),N);
for j=1:size(C.dat,2)
  if mod(j,10)==0
    disp(j);
  end
  
  % make sure that b is not above the whole range of the data
  r=range(C.dat(:,j));
%  mu=mean(C.dat(:,j));
%  sig=std(C.dat(:,j));
%  wgt=10./(10+abs((C.dat(:,j)-mu)/sig).^4);
%  wgt=ones(size(C.dat,1),1);
%  wgt=wgt./sum(wgt);
  if (1)
    for bi=1:N
      b=rng(1)+(bi-1)*stp;
      if b<=r
        xx=C.dat(:,j);
        xx=repmat(xx,1,10)+normrnd(0,0.05,size(xx,1),10);
        pos=mod(xx(:),b)/b;
        s(j,bi)=var(pos);
%        pos=mod(C.dat(:,j),b)/b;
%        s(j,bi)=sum(wgt.*pos.^2)-(sum(wgt.*pos))^2;
      else
        break;
      end
    end
    [ms,msi]=min(s(j,:));
    M(j)=ms;
    B(j)=(msi-1)*stp+rng(1);
  else
    if r<=rng(1)
      B(j)=r;
      M(j)=0;
    else
      [b,ms,status]=fminbnd(@(x) var(mod(C.dat(:,j),x)/x),rng(1),min(rng(2),r));
      B(j)=b;
      M(j)=ms;
    end
  end
  A(j)=mean(mod(C.dat(:,j),B(j)));
end
Z=C;
Z.dat=(C.dat-2)./repmat(B,size(C.dat,1),1)+2;

