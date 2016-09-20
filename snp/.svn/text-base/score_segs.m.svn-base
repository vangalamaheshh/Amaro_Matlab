function scores = score_segs(Q,alpha,lambda)
  
  %% Assumes probability function follows an exponential given by Pr(a) =
  %alpha*exp(-alpha*a)
  
  % If so, then score = -log(Pr(a)) = alpha*a - log(alpha)-log(lambda) =
  % alpha*a - log(alpha*lambda)
    
  scores = alpha*Q(:,12)-log(alpha*lambda);
  
  
