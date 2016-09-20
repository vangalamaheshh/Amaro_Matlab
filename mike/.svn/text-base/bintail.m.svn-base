function p = bintail(n,N,mu)
% returns value of 1-binocdf(n-1,N,mu)
% using incomplete beta function for increased precision

  if any(n<0), error('bintail: please supply n *not* n-1'); end
  N = max(n,N);

%  p = 1-binocdf(n-1,N,mu);

  p = betainc(1-mu,N-n+1,n,'upper');

  p(n==0) = 1;
  p(N==0) = 1;
