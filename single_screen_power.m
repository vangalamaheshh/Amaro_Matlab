function [po,ks]=single_screen_power(r0,r1,n,Ng,rnk,q,L)
if ~exist('q','var')
  q=0.1;
end
%r0=1.2e-6*1500*3;
%r1=r0+0.1;

%n1=11;
%n2=24;
%Ng=13023;
%rnk=6;
p=single_screen_prob(0:20,n,r0,'binom',L);
if p(end)>1e-8
  keyboard
end
pv=1-cumsum([0 p]);
ks=min(find(pv<=q/Ng*rnk));
ks=ks-1;

p1=single_screen_prob(0:20,n,r1,'binom',L);
po=1-sum(p1(1:ks)); 
