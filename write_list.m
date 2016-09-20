function write_list(fname,lst)

f=fopen(fname,'w');
if ischar(lst)
  lst=cellstr(lst);
end
for i=1:length(lst)
  fprintf(f,'%s\r\n',lst{i});
end

fclose(f);
