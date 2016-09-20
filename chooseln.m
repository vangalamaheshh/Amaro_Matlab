function res=chooseln(x,n)

if (n<0) | (x<0) | (n-x < 0)
    error('illegal numbers to chooseln');
end

res=gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1);

