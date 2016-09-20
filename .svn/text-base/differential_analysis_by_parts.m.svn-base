function [P,S]=differential_analysis_by_parts(D,cls0,cls1,test_type,nparts,lsfdir)

if nargin<5
  nparts==1;
end
if nargin<6
  lsfdir='./';
end

n=size(D.dat,1);
p=get_parts(1:n,nparts);

if nparts>1
  l=lsf(lsfdir);
  h=zeros(nparts,1);
  for i=1:nparts
    Di=[];
    Di.dat=D.dat(p{i},:);
    [l,h(i)]=bsub(l,{'P','S'},'differential_analysis',{Di,cls0,cls1,test_type});
  end
  [l,res]=wait(l); % wait for all
  P=NaN*ones(n,1);
  S=NaN*ones(n,1);
  for i=1:nparts
    P(p{i})=res{h(i)}.P;
    S(p{i})=res{h(i)}.S;    
  end
else
  [P,S]=differential_analysis(D,cls0,cls1,test_type);
end



