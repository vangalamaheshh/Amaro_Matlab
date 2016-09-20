function st=make_sup_list(name,strs)

st=[ name ': '];
if ischar(strs)
  strs=cellstr(strs);
end

for i=1:length(strs)
  st=[st num2str(i) '-' strs{i}];
  if (i<length(strs))
    st=[st '/'];
  end
end

