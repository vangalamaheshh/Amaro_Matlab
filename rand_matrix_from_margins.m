function M = rand_matrix_from_margins(A,nperm)
  % Generates uniform m x n random matrices with same marginal row and
  % column sums as the m x n matrix A
    
  % M = rand_matrix_from_marings(R,C)  
  % Inputs: A = m x n matrix input
  %         nperm (default 1) = number of random matrices to return
  %  
  % Outputs: M = m x n x nperm random matrix; each 2-D m x n matrix
  % indexed by first two dimensions has marginal row sums R and column
  % sums C
    
    if ~exist('nperm','var') || isempty(nperm)
      nperm = 1;
    end
    
    m = size(A,1);
    n = size(A,2);
    
    R = sum(A,2);
    C = sum(A,1);
    
    [mx mx_row] = max(R);
    min_rows = setdiff(1:m,mx_row);
    
    N = sum(R);
    N1 = N-R(mx_row);
        
    M = zeros(m,n,nperm);
            
    for j=1:nperm
      X = zeros(2,sum(R));
      [row1 filled] = rowperm(N,R,mx_row,min_rows);
            
      X(1,:) = row1;
      
      row2 = [];
      for i=1:n
        row2 = [row2 repmat(i,1,C(i))];
      end
      
      X(2,:) = row2;
      
      yy = zeros(m,n);
      yy(min_rows,1:max(X(2,filled))) = crosstab2d(X(1,filled)',X(2,filled)');
      yy(mx_row,:) = C-sum(yy,1);
      
      M(:,:,j) = yy;
      
    end
    
