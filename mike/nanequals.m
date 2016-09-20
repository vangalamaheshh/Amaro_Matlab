function eq = nanequals(a,b)
eq = (a==b) | (isnan(a)&isnan(b));
