function prob = mv_hypergeometric(A,R,C)
  % Computes the multivariate hypergeometric point probability of
  % observing a set of values (Aij) conditional on row sum R and column
  % sum C
    
    
  [m,n] = size(A);
  if length(R) ~= m
    error('R must have same number of elements as rows of A')
  end
  
  if length(C) ~= n
    error('C must have same number of elements as columns of A')
  end
  
  if sum(C) ~= sum(R)
    error('C and R should have same total sums!')
  else
    N = sum(R);
  end
  
  prob = prod(factorial(R))*prod(factorial(C))/(factorial(N)* ...
                                                prod(prod(factorial(A))));
  
