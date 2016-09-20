function remove_cols_rows(rows,cols,fname,fout)

fid=fopen(fname,'r');
fo=fopen(fout,'w');
ln=1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if(~ismember(ln,rows))
      colnum=cumsum(tline==9)+1; % find the tabs
      fprintf(fo,'%s\n',tline(find(~ismember(colnum,cols))));
  end
  ln=ln+1;
end
fclose(fo);
fclose(fid);


