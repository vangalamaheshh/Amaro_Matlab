function r=replace_by_rank(mat)

r=zeros(size(mat,1),size(mat,2));
for i=1:size(mat,2)
  [sv,si]=sort(mat(:,i));
  [ssi,ri]=sort(si);
  r(:,i)=ri;
end
