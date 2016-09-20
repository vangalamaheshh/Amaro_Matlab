function n=line_count(fname)

n=0;
if ispc
    sz=5000000;
else
    sz=250000000;
end
f=fopen(fname,'r');
a=zeros(sz,1);
if f>0
  while ~feof(f)
    [a,count]=fread(f,sz,'uchar');
    n=n+length(find(a==10));
  end
  verbose(['... ' num2str(n) ' lines read (line_count.m)'],30)
else
  error(['File not found: ' fname]);
end
fclose(f);
