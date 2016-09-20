function c=strs2cell(strs)

c=cell(size(strs,1),1);
for i=1:size(strs,1)
  c{i}=deblank(strs(i,:));
end
