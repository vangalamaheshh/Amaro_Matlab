function [g,val,members]=lnk_to_graph(lnk,N)
% assume lnk goes from leaves to root

br=1:N;
dad=zeros(2*N,1);
val=zeros(2*N,1);
members=zeros(2*N,2);
members(:,1)=N+1;
members(:,2)=-1;
members(1:N,1)=(1:N)';
members(1:N,2)=(1:N)';
dad(1:N)=1:N;
val(1:N)=NaN;
cur_node=N;
last=NaN;
for i=1:size(lnk,1)
  cur=lnk(i,:);
  node1=br(cur(1));
  node2=br(cur(3));
  if cur(5) ~= val(node1) & cur(5) ~= val(node2)
    cur_node = cur_node +1;
  end
  members(cur_node,:)=[min(members(cur_node,1),cur(1)) max(members(cur_node,2),cur(4))];
  dad(node1)=cur_node;
  dad(node2)=cur_node;
  br(cur(1))=cur_node;
  val(cur_node)=cur(5);
end
dad=dad(1:(cur_node-1));
g=sparse(1:(cur_node-1),dad,1:(cur_node-1),cur_node,cur_node);
members=members(1:cur_node,:);
val=val(1:cur_node);
