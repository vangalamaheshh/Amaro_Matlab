function d=binoci_diff(ci,x,n,alpha)

d(1)=abs(diff(betapdf(ci,x+1,n-x+1)));
d(2)=abs(diff(betacdf(ci,x+1,n-x+1)))-(1-alpha);
