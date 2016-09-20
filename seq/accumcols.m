function res=accumcols(subs,vals)

tmp=accumarray(subs,vals(:,1));
res=zeros(size(tmp,1),size(vals,2));
res(:,1)=tmp;
for i=2:size(vals,2)
	res(:,i)=accumarray(subs,vals(:,i));
end


