function v=frequent_vals(x,n)

[uu,ui,uj]=unique(x);
h=histc(uj,1:length(uu));
v=uu(find(h>=n));
