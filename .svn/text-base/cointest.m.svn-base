function p=cointest(n,N,f)

if (n>0)
  pL=binocdf(n-1,N,f);
else
  pL=0;
end

pE=binopdf(n,N,f);
pG=1-binocdf(n,N,f);

if pL<pG
  p=pE+pL;
  x=binoinv(p,N,1-f);
  if binocdf(x,N,1-f)==p
    verbose(num2str([ p binocdf(x,N,1-f) ]),1);
    p=2*p;
  else % it is greater
    if x>0
      verbose(num2str([p p+binocdf(x-1,N,1-f) p+binocdf(x,N,1-f)]),1);
      p=p+binocdf(x-1,N,1-f);
    end
  end
else
  p=pE+pG;
  x=binoinv(p,N,f);
  if binocdf(x,N,f)==p
    verbose([ p binocdf(x,N,1-f) ],1);
    p=2*p;
  else % it is greater
    if x>0
      verbose(num2str([p p+binocdf(x-1,N,f) p+binocdf(x,N,f)]),1);
      p=p+binocdf(x-1,N,f);
    end
  end
end
