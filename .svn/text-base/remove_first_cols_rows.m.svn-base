function remove_first_cols_rows(nrows,ncols,fname,fout)

fid=fopen(fname,'r');
fo=fopen(fout,'w');
ln=1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if(ln>nrows)
      pos=find(tline==9); % find the tabs
      fprintf(fo,'%s\n',tline((pos(ncols)+1):end));
  end
  ln=ln+1;
end
fclose(fo);
fclose(fid);


