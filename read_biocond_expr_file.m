function D=read_biocond_expr_file(fname)

gacc={};
sdesc={};

fid=fopen(fname,'r');
%F=fread(fid);
%s=char(F');
ln=1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  l{ln}=tline;
  ln=ln+1;
  if mod(ln,1000)==0
    disp(ln)
  end
end
fclose(fid);
ln=ln-1;

if ln<1
  error('FILE must have at least 1 lines');
end

sd=dlmsep(l{1});

sdesc=[];
for i=1:length(sd)
  sdesc=strvcat(sdesc,sd{i});
end

ngenes=size(l,2)-1;
dat=zeros(ngenes,size(sdesc,1));

gacc={};
for i=2:(ngenes+1)
  if mod(i,1000)==0
    disp(i)
  end
  la=dlmsep(l{i});
  gacc{i-1}=la{1};
  dat(i-1,:)=str2num(str2mat(la(2:end)))';
end

D.gacc=gacc;
D.sdesc=sdesc;
D.dat=dat;


