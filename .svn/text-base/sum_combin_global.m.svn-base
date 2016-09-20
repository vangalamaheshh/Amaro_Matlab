function y=sum_combin_global(n,x)
% how many ways to break n into x bins

global SUM_COMBIN

if isempty(SUM_COMBIN) || length(SUM_COMBIN)<x || length(SUM_COMBIN{x})<(n+1)
  disp(['Calculating ' num2str(n) ' ' num2str(x) ]);
  y=sum_combin(n,x);
  SUM_COMBIN{x}{n+1}=y;
else
  y=SUM_COMBIN{x}{n+1};
end
