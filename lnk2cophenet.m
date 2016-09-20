function m=lnk2cophenet(lnk,use_level)

if nargin==1
  use_level=0;
end

if use_level
  lnk(:,5)=flipud((1:size(lnk,1))');
end

N=max(max(lnk(:,1:4)));
m=zeros(N,N);

for i=1:size(lnk,1)
%  x=m(lnk(i,1):lnk(i,2),lnk(i,3):lnk(i,4)); 
  m(lnk(i,1):lnk(i,2),lnk(i,3):lnk(i,4))=lnk(i,5);
end
m=m+m';
