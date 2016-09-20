function y=binoconvpdf(x,n,p)

m=max(x);
for i=1:length(n)
  d{i}=binopdf(0:m,n(i),p(i));
end
dd=d{1};
for i=2:length(n)
  dd=conv(dd,d{i});
end
y=dd(x+1);
