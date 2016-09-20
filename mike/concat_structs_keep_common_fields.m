function X = concat_structs_keep_common_fields(X)

flds = fieldnames(X{1});
for i=2:length(X)
  ff = fieldnames(X{i});
  keep = true(length(flds),1);
  for j=1:length(flds), if ~ismember(flds{j},ff), keep(j)=false; end, end
  flds = flds(keep);
end
for i=1:length(X)
  X{i} = keep_fields(X{i},flds);
end
X = concat_structs(X);
