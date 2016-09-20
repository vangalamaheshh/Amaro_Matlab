function [po,ks]=single_screen_power_smooth(r0,r1,n,Ng,rnk,q,L,MAX,use_gaussian)
if ~exist('q','var')
  q=0.1;
end
if ~exist('use_gaussian','var')
  use_gaussian=0;
end
%r0=1.2e-6*1500*3;
%r1=r0+0.1;

%n1=11;
%n2=24;
%Ng=13023;
%rnk=6;
if ~exist('MAX','var')
  MAX=50;
end

if ~use_gaussian
  p=single_screen_prob(0:MAX,n,r0,'binom',L);
else
  fprintf(1,'.');
  p=single_screen_prob(0:MAX,n,r0,'gaussian',L);
end
if p(end)>1e-8
  keyboard
end

pv=1-cumsum([0 p]);
% pv1=[1 1-binocdf(0:MAX,n*L,r0/L)]; % OK up to 1e-8

sig=q/Ng*rnk;
ks=min(find(pv<=sig));
if ks>1
  c=(sig-pv(ks))/(pv(ks-1)-pv(ks));
else
  c=0;
end
ks=ks-1;

if ~use_gaussian
  p1=single_screen_prob(0:MAX,n,r1,'binom',L);
else
  fprintf(1,'.');
  p1=single_screen_prob(0:MAX,n,r1,'gaussian',L);
end
po=1-sum(p1(1:ks)); 
po=po+c*p1(ks);



