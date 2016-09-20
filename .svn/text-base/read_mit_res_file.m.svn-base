function D=read_mit_res_file(fname)

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

if ln<4
  error('RES FILE must have at least 3 lines');
end

sd=dlmsep(l{1});
% Description \t Accession Name1 \t \t Name2 \t \t Name\3 \t

if ~strcmp(sd{1},'Description')
  error(['RES FILE first element of first row should be ' ...
         'Description']);
end
if ~strcmp(sd{2},'Accession')
  error(['RES FILE second element of first row should be ' ...
         'Accession']);
end

sdesc=[];
for i=3:2:length(sd)
  sdesc=strvcat(sdesc,sd{i});
end

la=dlmsep(l{2});
sscale=[];
for i=3:2:length(la);
  sscale=strvcat(sscale,la{i});
end

ngenes=sscanf(l{3},'%d');
if ngenes+3 ~= length(l)
  error('RES FILE inconsistent number of lines and genes');
end

dat=zeros(ngenes,size(sdesc,1));
affy_call=zeros(ngenes,size(sdesc,1));

gdesc={};
gacc={};

for i=4:(ngenes+3)
  if mod(i,1000)==0
    disp(i)
  end
  la=dlmsep(l{i});
  gdesc{i-3}=la{1};
  gacc{i-3}=la{2};
  dat(i-3,:)=str2num(str2mat(la(3:2:end)))';
  temp=str2mat(la(4:2:end));
  temp(temp=='A')=0;
  temp(temp=='P')=2;
  temp(temp=='M')=1;
  temp=double(temp);
  affy_call(i-3,:)=temp';
end

D.gdesc=gdesc;
D.gacc=gacc;
D.sdesc=sdesc;
D.sscale=sscale;
D.dat=dat;
D.affy_call=affy_call;
