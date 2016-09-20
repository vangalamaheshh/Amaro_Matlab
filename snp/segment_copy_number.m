function [seg,p,t]=segment_copy_number(dat,sz,diff_method,pcutoff)

hsz1=floor(sz/2);
hsz2=sz-hsz1;

p=zeros(1,length(dat)-1); 
t=zeros(1,length(dat)-1); 
seg=zeros(1,length(dat)-1); % break after point i

old=0;

if old
  for i=1:(length(dat)-1)
    g1=dat(max(i-hsz1+1,1):i);
    g2=dat((i+1):min(i+hsz2,length(dat)));
    D.dat=[g1' g2'];
    [p(i),t(i)]=differential_analysis(D,1:length(g1),length(g1)+(1:length(g2)),diff_method,0);
  end

  seg=double(p<=pcutoff);
  seg=cumsum([1 seg]);
  
  return
end

for i=1:min(hsz1,length(dat)-1)
  g1=dat(max(i-hsz1+1,1):i);
  g2=dat((i+1):min(i+hsz2,length(dat)));
  D.dat=[g1' g2'];
  [p(i),t(i)]=differential_analysis(D,1:length(g1),length(g1)+(1:length(g2)),diff_method,0);
end

% keyboard
x=zeros(hsz1+hsz2,length(dat));
for j=1:(hsz1+hsz2)
  x(j,1:(length(dat)-j+1))=dat(j:end);
end
D.dat=x(:,2:(length(dat)-hsz1-hsz2))';
cls0=1:hsz1;
cls1=hsz1+(1:hsz2);
[p((hsz1+1):(length(dat)-hsz2-1)),t((hsz1+1):(length(dat)-hsz2-1))]=differential_analysis(D,cls0,cls1,diff_method,0);

for i=max((length(dat)-hsz2),1):(length(dat)-1)
  g1=dat(max(i-hsz1+1,1):i);
  g2=dat((i+1):min(i+hsz2,length(dat)));
  D.dat=[g1' g2'];
  [p(i),t(i)]=differential_analysis(D,1:length(g1),length(g1)+(1:length(g2)),diff_method,0);
end


segup=double(p<=pcutoff & t>0);
segdown=double(p<=pcutoff & t<0);

segup_pos=find(segup);
[ss,si]=sort(p(segup_pos));
segup_pos=segup_pos(si);
for i=1:(length(segup_pos)-1)
  segup_pos(i+find(abs(segup_pos((i+1):end)-segup_pos(i))<hsz1))=-1;
end
segup_pos(segup_pos==-1)=[];

segdown_pos=find(segdown);
[ss,si]=sort(p(segdown_pos));
segdown_pos=segdown_pos(si);
for i=1:(length(segdown_pos)-1)
  segdown_pos(i+find(abs(segdown_pos((i+1):end)-segdown_pos(i))<hsz1))=-1;
end
segdown_pos(segdown_pos==-1)=[];

seg=zeros(1,length(p));
seg(segup_pos)=1;
seg(segdown_pos)=1;

seg=cumsum([1 seg]);

