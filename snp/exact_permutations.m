function [thst_amp,thst_del]=exact_permutations(C,t1s,t2s)


xamp=C.smooth;
xamp(xamp<t2s)=0;

for i=1:size(xamp,2)
    ha{i}=histc(xamp(:,i),0:0.01:10)/size(xamp,1);
end 

camp=ha{1};
fprintf(1,'Amp: 1 ',i);
for i=2:length(ha);
    camp=conv(camp,ha{i});
    fprintf(1,'%d ',i);
end
fprintf(1,'\n');
if length(camp)<10010
  tmp=zeros(10010,1);
  tmp(1:length(camp))=camp;
  camp=tmp;
end

thst_amp=sum(reshape(camp(1:10010),10,1001),1)';

xdel=-C.smooth;
xdel(xdel<-t1s)=0;

for i=1:size(xdel,2)
    hd{i}=histc(xdel(:,i),0:0.01:10)/size(xdel,1);
end 

cdel=hd{1};
fprintf(1,'Del: 1 ',i);
for i=2:length(hd);
    cdel=conv(cdel,hd{i});
    fprintf(1,'%d ',i);
end
fprintf(1,'\n');

if length(cdel)<10010
  tmp=zeros(10010,1);
  tmp(1:length(cdel))=cdel;
  cdel=tmp;
end
thst_del=sum(reshape(cdel(1:10010),10,1001),1)';

