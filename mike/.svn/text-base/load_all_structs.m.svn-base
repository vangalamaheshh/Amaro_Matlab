function x = load_all_structs(filemask)

fprintf('Getting list of %s\n',filemask);
d = direc(filemask);
fprintf('%d files found.\n',length(d));
x = cell(length(d),1);
for i=1:length(d)
  fprintf('Loading %d/%d    \t%s\n',i,length(d),d{i});
  x{i} = load_struct(d{i});
  x{i}.sourcefile = repmat({d{i}},slength(x{i}),1);
end

