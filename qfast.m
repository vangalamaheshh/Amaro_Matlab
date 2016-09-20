function Q=qfast(p,n)

nz=find(p>0);
q=p(nz);
if n>1
  x=-log(q);
end
t=q;
for i=1:(n-1)
  t=t.*x/i;
  q=q+t;
end
Q=p;
Q(nz)=q;

