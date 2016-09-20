function [D,A,probesets]=read_novartis_data(fname,cells)

dlm=',';
% Probe Set Name,ID,Gene,cellname,pname,panelnbr,cellnbr,Signal,Detection,P Value
% 36460_at,GC26855_A,POLR1C,K-562,Leukemia,7,5,95.758659,P,0.02786
if exist([ fname '_small.mat'],'file')
  load([ fname '_small.mat']);
else
  fid=fopen([ fname '.txt' ]);
  C=textscan(fid,'%s%*s%*s%*s%*s%d%d%f%s%f','delimiter', ',','headerLines',1);
  fclose(fid);
  save([ fname '_small.mat'],'C');
end

clt=double(C{2})+double(C{3})/100;
cells2=cells(:,1)+cells(:,2)/100;

idx=find(ismember(clt,cells2));
clt2=clt;
for i=1:length(cells2)
  clt2(clt==cells2(i))=i;
end

[probesets,probei,probej]=unique(strvcat(C{1}),'rows');
nps=size(probesets,1);
I=ones(nps,length(cells));
D=zeros(nps,length(cells),3);
A=zeros(nps,length(cells),3);
return

for i=1:length(idx)
  ci=idx(i);
  clid=clt2(ci);
  pr=probej(ci);
  D(pr,clid,I(pr,clid))=C{4}(ci); 
  A(pr,clid,I(pr,clid))=convert_enum(C{5}{ci},{'A',-1;'P',1;'M',0});
  I(pr,clid)=I(pr,clid)+1;
  if mod(i,10000)==0
    disp(i);
  end
end

