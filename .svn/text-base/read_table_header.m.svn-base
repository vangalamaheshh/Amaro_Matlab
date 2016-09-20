function D=read_table_header(fname,dlm,nhead,nleft,ignore_lines)
f=fopen(fname,'r');
for i=1:ignore_lines
  first_line=fgetl(f);
end
first_line=fgetl(f);
ncol=length(find(first_line==dlm))+1;
fclose(f);

f=fopen(fname,'r');
r=textscan(f,repmat('%s',1,ncol),nhead,'delimiter',dlm,'bufSize',1000000,'headerLines',ignore_lines);
r2=textscan(f,[repmat('%s',1,nleft) repmat('%n',1,ncol-nleft)],'treatAsEmpty','na','delimiter',dlm,'bufSize',1000000);
fclose(f);

D.sdesc=cat(1,r{(1+nleft):end});
D.gacc=r2{1};
D.dat=cat(2,r2{(1+nleft):end});
D.gdesc=r2{2};

