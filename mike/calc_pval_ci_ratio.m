function ci_ratio = calc_pval_ci_ratio(k,n)
% ci_ratio = calc_pval_ci_ratio(k,n)
%
% k = # successes
% n = # trials
%
% based on manuscript by Gaddy Getz

k = max(0,k);

ci_ratio = 2*1.96*sqrt((n-k+1) ./ ((k+1) .* (n+3)));

