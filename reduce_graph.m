function redg=reduce_graph(g)
% x-->y-->z1  ==> -->y-->z1
%     |-->z2         |-->z2
% assume only one dad to each node

oneson=find(sum(spones(g),1)==1);
for i=1:length(oneson)
  [son,dum,dum2]=find(g(:,oneson(i))); % find father of oneson(i)
  g(son,oneson(i))=0; 
  [dum,dad,dum2]=find(g(oneson(i),:));
  g(son,dad)=1;
  g(oneson(i),:)=0;
end

redg=g;
