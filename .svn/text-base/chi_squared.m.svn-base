function [chis ps flags] = chi_squared(A)
  % Calculates the pearson chi-squared value (and p-value) for m x n x p matrix A
  % Returns chis (p x 1 vector for each m x n submatrix of A)
  %         ps (p x 1 vector of p-values for each chi)
  %         flags (p x 1 vector of warning flags for each matrix).  
  %         Flag = 1 if matrix is very inhomogeneous [abs(E-A) > E for
  %         some matrix element] and g-test is preferred
  %              = 2 if any E is less than 10 (and exact or Monte Carlo
  %              method is preferred).
    
  flags = zeros(1,size(A,3));
  
  for j=1:size(A,3)
    R = sum(A(:,:,j),2);
    C = sum(A(:,:,j),1);
    N = sum(R);
    
    E = R*C/N;
  
    chis(j) = sum(sum((E-A(:,:,j)).^2./E));
  
    if any(any(abs(E-A(:,:,j)) > E))
      flags(j) = 1;
      verbose('Chi-squared p-value not accurate when abs(Ei-Ai) > Ei for any bin.  Use g-test instead.',30);
    end
    
    if any(any(E)) < 10
      flags(j) = 2;
      verbose('Chi-squared p-value not accurate when any expected bin < 10!',30)
    end
  
    nu = (size(A(:,:,j),1)-1)*(size(A(:,:,j),2)-1);
  
    ps(j) = 1 - chi2cdf(chis(j),nu);
  end
  
