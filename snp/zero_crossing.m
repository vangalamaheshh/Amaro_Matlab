function n=zero_crossing(C)

pos=C.dat>0;
neg=C.dat<0;
pn=sum(pos(1:(end-1),:).*neg(2:end,:));
np=sum(neg(1:(end-1),:).*pos(2:end,:));
n=pn+np;
