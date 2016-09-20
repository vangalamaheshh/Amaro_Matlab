function a=make_cell_array(c)

w=length(c{1});
a=cell(length(c),w);
for i=1:length(c)
  for j=1:(min(length(c{i}),w))
    a{i,j}=c{i}{j};
  end
end
