function mst=mstree(sedges,n)

lb=1:n;
mst=zeros(n-1,3);
j=1;
for i=1:size(sedges,1)
	if j == n 
		break
	end

	if  lb(sedges(i,1)) ~= lb(sedges(i,2)) 
		mst(j,:)=sedges(i,:);
		a=find(lb==lb(sedges(i,2)));
		lb(a)=lb(sedges(i,1))*ones(length(a),1);
		j = j+1;
	end
end

mst=mst(1:(j-1),:);
