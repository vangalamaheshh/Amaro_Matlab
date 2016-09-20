function [class,score,knnidx]=knn_classify(pts,mempts,cls,w,k,dist_type)
k=min(k,size(mempts,1));
knnidx=zeros(k,size(pts,1));
d=dist(mempts,pts,dist_type);
for i=1:size(pts,1)
%  d=dist(mempts,pts(i,:),dist_type);
%  [ds,di]=sort(d);
  [ds,di]=sort(d(:,i));
  knnidx(:,i)=di(1:k);
  ck=sum(cls(:,di(1:k)),2).*w;
  [score(i),class(i)]=max(ck);
end

