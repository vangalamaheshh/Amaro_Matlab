function z=lnk2z(lnk,idx)
% z=lnk2z(lnk,idx)
%    convert lnk,idx to a z matrix (output of linkage)
%    z is in the format  
%         node1 node2 dist12
%         node3 node4 dist34
%    lnk format is: to from to from dist
%    idx is the reorder indices
% 

n=length(idx);
z=zeros(size(lnk,1),3);
nd=idx;
for i=1:size(lnk,1)
  set1=lnk(i,1):lnk(i,2);
  set2=lnk(i,3):lnk(i,4);
  z(i,:)=[min(nd(set1)) min(nd(set2)) lnk(i,5)];
  nd([set1 set2])=i+n;
end

