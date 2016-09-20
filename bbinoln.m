function lp=bbinoln(k,n,a,b)

    lp =gammaln(a+b)-gammaln(a)-gammaln(b) + ...
         gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
         gammaln(a+k) + gammaln(b+n-k) - gammaln(n+a+b);         


