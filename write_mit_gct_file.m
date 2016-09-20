function write_mit_gct_file(fname,D,format_str)
if ~exist('format_str','var')
  format_str=[];
end
tab=char(9);
fid=fopen(fname,'w');
fprintf(fid,'#1.2\n');
ngenes=size(D.dat,1);
nsamples=size(D.dat,2);

fprintf(fid,'%d\t%d\n',ngenes,nsamples);
    fprintf(fid,'%s%s%s%s%s\n','Name',tab,'Description',tab,write_line(D.sdesc,tab));

for i=1:ngenes
  if mod(i,1000)==0
    disp(i);
  end
  if iscell(D.gacc)
    fprintf(fid,'%s%s',D.gacc{i},tab);
  else
    fprintf(fid,'%s%s',deblank(D.gacc(i,:)),tab);
  end
  if iscell(D.gdesc)
    fprintf(fid,'%s%s',D.gdesc{i},tab);
  else
    fprintf(fid,'%s%s',deblank(D.gdesc(i,:)),tab);
  end
%  x=write_line(num2str(D.dat(i,:)'),tab);
%  fprintf(fid,'%s\n',x);
  fprintf(fid,[ repmat(['%' format_str 'f\t'],1,size(D.dat,2)-1) '%' format_str 'f\n'],D.dat(i,:));
end

fclose(fid);



