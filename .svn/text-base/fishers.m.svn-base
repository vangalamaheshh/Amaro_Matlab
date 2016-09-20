function p = fishers(M)
  %% Performs fisher's exact test on the m x n contingency table
  %represented by M
    
    [m,n] = size(M);
    
    R = sum(M,2); % row sums
    C = sum(M,1); % col sums
    
    N = sum(R); % total
           
    As = enum_mat_with_fixed_margins(R,C);
    
    % Find index corresponding to M
    
    F = @(xx) prod(prod(double(eq(xx,M))));
    idx = find(cellfun(F,As));

    prob = [];
    chi = [];
    for j=1:length(As)
      chi(j) = chi_squared(As{j});
      prob(j) = mv_hypergeometric(As{j},R,C);
    end
    
    [C csi] = sort(chi,'ascend');
    c_idx = find(csi == idx);
    
    p = sum(prob(csi(c_idx:end)));
    
