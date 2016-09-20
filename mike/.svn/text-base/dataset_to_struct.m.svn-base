function s = dataset_to_struct(d)

p = get(d,'VarNames');
s = [];
for i=1:length(p)
  s = setfield(s,p{i},double(d(:,i)));
end
