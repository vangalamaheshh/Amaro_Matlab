function eq = naneq(a,b)
eq = (a==b) | (isnan(a)&isnan(b));
