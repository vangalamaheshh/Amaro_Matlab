function y=ln_binocdf(x,n,p,tail)

if ~exist('tail','var')
  tail=-1;
end

f=ln_binopdf(0:n,n,p);
if prod(size(x))>(n+1)
  lcdf=ln_binocdf(0:n,n,p);
  y=lcdf(x+1);
else
  y=zeros(size(x));
  for i=1:prod(size(x))
    if tail==-1 % calc left tail
      if x(i)==0
        y(i)=f(1);
      elseif x(i)==n
        y(i)=0;
      else
        r=1:(x(i)+1);
        [mf,mfi]=max(f(r));
        y(i)=mf+log(sum(exp(f(r)-mf)));
      end
    else
      if x(i)==0
        y(i)=0;
      elseif x(i)==n
        y(i)=f(n+1);
      else
        r=(x(i)+1):(n+1);
        [mf,mfi]=max(f(r));
        y(i)=mf+log(sum(exp(f(r)-mf)));
      end
    end
  end
end

% [ log(binocdf(0:100,100,0.01)); ln_binocdf(0:100,100,0.01) ]
% [ log([ 1 1-binocdf(0:99,100,0.01)]); ln_binocdf(0:100,100,0.01,1) ] % right tail
