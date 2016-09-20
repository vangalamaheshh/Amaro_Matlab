function C = char2cell(c)

C = cell(size(c,1),1);
for i=1:size(c,1)
  C{i} = c(i,:);
end
