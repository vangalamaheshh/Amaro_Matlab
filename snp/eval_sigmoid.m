function F = eval_sigmoid(x,theta)
  
%% Evaluates the function F = vmax*x^n/(k^n+x^n)
%  where theta = [vmax k n]
%  x can be a scalar, vector or matrix
    
  vmax = theta(1);
  k = theta(2);
  n = theta(3);
  
  F = (vmax*x.^n)./(k^n+x.^n);
