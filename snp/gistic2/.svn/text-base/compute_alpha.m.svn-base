function [alpha beta] = compute_alpha(hc,hx)
  
  [A i1] = min(abs(hx - quantile(hx,.25)));
  [A i2] = min(abs(hx - quantile(hx,.75)));

  p = polyfit(hx(i1:i2)',log(hc(i1:i2)),1);
  
  alpha = -1*p(1);
  
  beta = exp(p(2));
     
  
