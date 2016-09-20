function [low,high]=two_sided_prctile(x,f)


[u,a,b]=unique(x);
n=floor(length(x)*f);
low=1;
high=length(u);
h=hist(b,1:length(u));
stuck=0;
t=0;
while (t <= n) & ~stuck
  if (h(high)==h(low)) & (t+h(high)+h(low)<=n)
    high=high-1;
    low=low+1;
    t=t+h(high)+h(low);
  elseif (h(high)<=h(low)) & (t+h(high)<=n)
    high=high-1;
    t=t+h(high);
  elseif t+h(low)<=n
    low=low+1;
    t=t+h(low);
  else
    stuck=1;
  end
end

low=u(low);
high=u(high);

%%%%% see fix_Pvalues2 for ideas of fast ways to doing so

