function D=read_mit_gct_file2(fname)

gdesc={};
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


if ln<3
  error('GCT FILE must have at least 2 lines');
end

% skip line 1,2

sd=dlmsep(l{3});
% Description \t Accession Name1 \t \t Name2 \t \t Name\3 \t

if ~strcmp(sd{1},'Name')
  error(['GCT FILE first element of the third row should be ' ...
         'Name']);
end
if ~strcmp(sd{2},'Description')
  error(['GCT FILE first element of the third row should be ' ...
         'Description']);
end

sdesc=[];
for i=3:length(sd)
  sdesc=strvcat(sdesc,sd{i});
end

tmp=sscanf(l{2},'%d\t%d');
ngenes=tmp(1);
nsamples=tmp(2);
if ngenes+3 ~= length(l)
  error('GCT FILE inconsistent number of lines and genes');
end

if nsamples ~= size(sdesc,1)
  error('GCT FILE inconsistent number of rows and samples');
end

dat=zeros(ngenes,size(sdesc,1));

gdesc={};
gacc={};

for i=4:(ngenes+3)
  if mod(i,1000)==0
    disp(i)
  end
  la=dlmsep(l{i});
  gdesc{i-3}=la{2};
  gacc{i-3}=la{1};
%  dat(i-3,:)=str2num(str2mat(la(3:1:end)))';
end

dat=dlmread(fname,'\t',4,3);

D.gdesc=gdesc;
D.gacc=gacc;
D.sdesc=sdesc;
D.dat=dat;


