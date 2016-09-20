function y=tpfpdf(x,n,tau)
% conv of two -log(truncated p-values)
% 
% s=1e-4;
% qqq=tpfpdf(s:s:50,2,0.05)*s;
% sum(qqq)+tpfpdf(0,2,0.05) % = 1
%
%tic;
%s=1e-4;
%qqq=tpfpdf(s:s:50,3,0.05)*s;
%sum(qqq)+tpfpdf(0,3,0.05) % = 1
%toc
%

if numel(x)>1
   y=nan(size(x));
   for i=1:numel(x)
       y(i)=tpfpdf(x(i),n,tau);
   end
   return
end
   
a=-log(tau);

switch n
    case 2
        if x<0
            y=0;
        elseif x==0
            y=(1-tau)^2;
        elseif x<a
            y=0;
        elseif x<2*a
            y=2*(1-tau)*exp(-x);
        else % x>= 2*a
            y=2*(1-tau)*exp(-x)+(x-2*a)*exp(-x);
        end
        
    case 3
        if x<0
            y=0;
        elseif x==0
            y=(1-tau)^3;
        elseif x<a
            y=0;
        elseif x<2*a
            y=3*(1-tau).^2*exp(-x);
        elseif x<3*a
            y=3*(1-tau).^2*exp(-x) + ...
              3*(1-tau)*(x-2*a)*exp(-x);
        else % x>= 3*a
            y=3*(1-tau).^2*exp(-x) + ...
              3*(1-tau)*(x-2*a)*exp(-x) + ...
              1/2*(x-3*a)^2*exp(-x);
        end
    otherwise
        error('not supported');
end

