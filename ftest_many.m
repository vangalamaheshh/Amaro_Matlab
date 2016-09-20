function [pv,Br,msr,dfr,Bf,msf,dff]=ftest_many(X,v1,v2)

v1=as_row(v1);
v2=as_row(v2);

% keyboard
pv=zeros(size(X,1),1);
if (0)
  for i=1:size(X,1)
    [p,t,stats,terms]=anovan(X(i,:),{v1,v2},'display','off','model','interaction');
    pv(i)=p(1);
    if mod(i,1000)==0
      disp(i);
    end
  end
end

% (SSR-SSF)/(dfR-dfF) / [ SSF/dfF ] 

% fit model 
% C coding of X
% B C = X => B = X inv(C);

%v1=(v1==1)-(v1==0);
%v2=(v2==1)-(v2==0);

N=length(v1);
%CF=[ 1 p c pc ];
% for 0/1 mapping use: CF=[ ones(N,1) v1' v2' (v1 & v2)'];
% 1/-1 mapping
CF=[ ones(N,1) v1'*2-1 v2'*2-1 ((v1*2-1).*(v2*2-1))'];
% for ttest use : CF=[ ones(N,1) v1'*2-1];
[ssf,Bf]=calc_SS(X,CF);
dff=N-size(CF,2);
msf=ssf./dff;

% CR=[ 1 c ];
% for 0/1 mapping use: CR=[ ones(N,1) v2' ];
CR=[ ones(N,1) v2'*2-1 ];
% for ttest use: CR=[ ones(N,1) ];
[ssr,Br]=calc_SS(X,CR);
dfr=N-size(CR,2);
delta_ssr=ssr-ssf;
msr=delta_ssr/(dfr-dff);
pv=1-fcdf(msr./msf,dfr-dff,dff);
pv=pv';

function [sse,B]=calc_SS(X,C)

B=pinv(C)*X';
Xe=C*B;
res=X'-Xe;
sse=sum(res.^2,1);
