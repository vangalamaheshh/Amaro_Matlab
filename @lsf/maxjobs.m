function res=maxjobs(l,m)

if nargin==1
  res=l.maxjobs;
else
  l.maxjobs=m;
  res=l;
end
