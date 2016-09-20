function  [P1sgte,P1slte,P2s]=fix_Pvalues_online(K,N,smooth_flag);

P1sgte=NaN*ones(size(N));
P1slte=P1sgte;
P2s=P1sgte;

if smooth_flag
  P1sgte=(K+1)./(N+2);
else
  P1sgte=(K+eps)./(N+eps);
end  
P1slte=1-P1sgte;


P2s=2*min(P1sgte,P1slte);

