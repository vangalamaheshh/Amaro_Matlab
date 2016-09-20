function write_mit_gmt_file(fname,gs)

f=fopen(fname,'w');
for i=1:length(gs)
  fprintf(f,'%s\t%s',gs(i).name,gs(i).type);
  for j=1:length(gs(i).genes)
    fprintf(f,'\t%s',gs(i).genes{j});
  end
  fprintf(f,'\n');
end
fclose(f);
