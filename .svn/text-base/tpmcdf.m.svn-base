function y=tpfcdf(x,n,tau)
% cdf of conv of two -log(truncated p-values)
% 

if numel(x)>1
   y=nan(size(x));
   for i=1:numel(x)
       y(i)=tpfcdf(x(i),n,tau);
   end
   return
end

a=-log(tau);
  
switch n
    case 2
        if x<0
            y=0;
        elseif x<a
            y=(1-tau)^2;
        elseif x<2*a
            y=(1-tau)^2 + ...
              2*tau*(1-tau)-2*(1-tau)*exp(-x);
        else % x>= 2*a
            y=(1-tau)^2 + ...
              2*tau*(1-tau)-2*(1-tau)*exp(-x) + ...
              tau.^2 - exp(-x) - x.*exp(-x) + 2*a*exp(-x);
        end
        
        
    case 3
        if x<0
            y=0;
        elseif x<a
            y=(1-tau)^3;
        elseif x<2*a
            y=(1-tau)^3 + ...
              3*(1-tau)^2*(tau-exp(-x));
        elseif x<3*a
            y=(1-tau)^3 + ...
              3*(1-tau)^2*(tau-exp(-x)) + ...
              3*exp(- 2*a - x)*(tau - 1)*(exp(2*a) - exp(x) - 2*a*exp(2*a) + x*exp(2*a));
        else % x>=3*a
            y= (1-tau)^3 + ...
              3*(1-tau)^2*(tau-exp(-x)) + ...
              3*exp(- 2*a - x)*(tau - 1)*(exp(2*a) - exp(x) - 2*a*exp(2*a) + x*exp(2*a)) + ...
              exp(-3*a) - (exp(-x)*(9*a^2 - 6*a*x - 6*a + x^2 + 2*x + 2))/2;
        end
        
        
    otherwise
        error('not supported');
end

