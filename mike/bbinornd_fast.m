function out = bbinornd_fast(N,a,b,mm,nn)

if ~exist('mm','var'), mm=1; end
if ~exist('nn','var'), nn=mm; end

%if (a<1 || b<1), error('requires (a,b) not (x,X)'); end

out = nan(mm,nn);

for mi=1:mm, for ni=1:nn
    
    r = rand;

    n = 0;
    z = bbinopdf(n,N,a,b);
    
    while z<r
      n=n+1;
      z=z+bbinopdf(n,N,a,b);
    end

    out(mi,ni)=n;
end,end


