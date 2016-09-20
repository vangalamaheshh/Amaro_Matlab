function X=collapse_identical_consecutive(X)

dat=X.dat;
dat(isnan(dat))=0;
dif=diff([-ones(1,size(dat,2)); dat],1);
st=find(any(dif,2));
en=[ st(2:end)-1; size(dat,1) ];
X.cdat=dat(st,:);
X.cdat(X.cdat==0)=NaN;

X.range=[st en];

