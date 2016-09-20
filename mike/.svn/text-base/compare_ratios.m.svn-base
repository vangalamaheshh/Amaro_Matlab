function p = compare_ratios(n1,N1,n2,N2)
% p = compare_ratios(n1,N1,n2,N2)

mu1 = n1/N1;
mu2 = n2/N2;
mustar = (n1+n2)/(N1+N2);

m = mustar*(1-mustar);
p = 1-normcdf(abs(mu1-mu2),0,sqrt((m/N1)+(m/N2)));


