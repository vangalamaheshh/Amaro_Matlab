function res=cellfun_any(fun,c)

res=cell(size(c,1),size(c,2));

if ~isempty(find(fun=='('))
  f=regexprep(fun,'x','c{i,j}');
else
  f=[ fun '(c{i,j})'];
end
for i=1:size(c,1)
  for j=1:size(c,2)
    eval(['res{i,j}=' f  ';']);
  end
  if mod(i,100)==0
    verbose(i);
  end
end

