function g = find_non_top_genes(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'include_silent_in_sort',false);
P=impose_default_value(P,'exclude_top_n_genes',40);

C = compute_mutation_rates(M,P);
g = C.nontop;
