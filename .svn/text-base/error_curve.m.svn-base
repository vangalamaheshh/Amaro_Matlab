function [err_c,xpos,bin]=error_curve(x,y,prc,nbins,min_per_bin)

xpos=min(x):(max(x)-min(x))/(nbins-1):max(x);
[n,bin]=histc(x,xpos);

for i=1:nbins
  idx=find(bin==i);
  if exist('min_per_bin','var') && length(idx)<min_per_bin
    err_c(i)=NaN;
  else
    err_c(i)=prctile(y(idx),prc);
  end
end

