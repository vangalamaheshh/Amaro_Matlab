function write_mit_res_file(fname,D)
tab=9;
fid=fopen(fname,'w');
ngenes=size(D.dat,1);
nsamples=size(D.dat,2);

fprintf(fid,'%s\r\n',[ 'Description' tab 'Accession' tab ...
        write_line(D.sdesc,[ tab tab ])  ]);
if isfield(D,'sscale') & ~isempty(D.sscale)
    fprintf(fid,'%s\r\n',[ tab ...
        write_line(D.sscale,[ tab tab ]) ]);
else
    fprintf(fid,'%s\r\n',[ tab ...
        write_line(repmat('NoDesc',size(D.dat,2),1),[ tab tab ])]);  
end

if ischar(D.gacc)
  D.gacc=cellstr(D.gacc);
end

if ischar(D.gdesc)
  D.gdesc=cellstr(D.gdesc);
end


%ngenes=100;
fprintf(fid,'%d\r\n',ngenes);
for i=1:ngenes
  if mod(i,1000)==0
    disp(i);
  end
  fprintf(fid,'%s',[D.gdesc{i} tab]);
  fprintf(fid,'%s',[D.gacc{i} tab]);
  
  %  comb=[ num2cell(D.dat(i,:))' num2affy_call(D.affy_call(i,:)')
  %  ]';
  
  if isfield(D,'affy_call')
    calls=num2affy_call(D.affy_call(i,:)');
  else
    calls=cellstr(repmat('M',size(D.dat,2),1));
  end
  for j=1:(size(D.dat,2)-1)
    fprintf(fid,'%f\t%s\t',D.dat(i,j),calls{j});
  end
  fprintf(fid,'%f\t%s\r\n',D.dat(i,end),calls{end});
end

fclose(fid);



