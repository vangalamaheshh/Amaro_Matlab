function [D,P,PA]=kol_smir(x1,x2,a,N)
if (size(x1,1)==1)
  x1=x1';
end
if (size(x2,1)==1)
  x2=x2';
end

if nargin>2
  x1=repmat(x1,1,N)+random('norm',0,a,size(x1,1),N);
  x1=x1(:);
  x2=repmat(x2,1,N)+random('norm',0,a,size(x2,1),N);
  x2=x2(:);
end

x=[x1; x2];
t=[ones(length(x1),1)/length(x1); -ones(length(x2),1)/length(x2)];

[sx,ix]=sort(x);
ct=cumsum(t(ix));

[D,Di]=max(abs(ct));

P=1;
PA=1;


