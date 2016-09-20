function res=discretize(vec,rngs)
res=NaN*ones(size(vec));
for i=1:length(rngs)
  res(find(vec>=rngs{i,1}(1) & vec < rngs{i,1}(2)))=rngs{i,2};
end

