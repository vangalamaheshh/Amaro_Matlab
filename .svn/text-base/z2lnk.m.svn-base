function [lnk,idx] = z2lnk( z )
% [lnk,idx] = z2lnk( z )
%    convert a z matrix (output of linkage) to lnk and idx
%    z is in the format  
%         node1 node2 dist12
%         node3 node4 dist34
%    lnk format is: to from to from dist
%    idx is the reorder indices
% 


m = size( z, 1 )+1;
idx=zeros(1,m);
j=2;
idx(1:2)=z(m-1,1:2);
while max( idx ) > m
  [k,i] = max( idx );
  idx(1:j+1) = [idx(1:i-1) z(k-m,1:2) idx(i+1:j)];
  j = j+1;
end
invidx=zeros(1,m);
invidx( idx ) = 1:m;

k=zeros(2,2);
lnk=zeros(m-1,5);
for i=1:m-1
  for j=1:2
    if z(i,j)<=m
      k(j,:)=[invidx(z(i,j)) invidx(z(i,j))];
    else
      k(j,:)=lnk( z(i,j)-m, [1 4] );
    end
  end
  if k(1,2)<k(2,1)
    lnk(i,:) = [k(1,:) k(2,:) z(i,3)];
  else
    lnk(i,:) = [k(2,:) k(1,:) z(i,3)];
  end
end

