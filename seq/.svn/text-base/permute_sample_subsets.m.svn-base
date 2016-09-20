function P_sig = permute_sample_subsets(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'permute_count',1000);
P=impose_default_value(P,'permute_size',24);
P=impose_default_value(P,'Q_threshold',0.1);

sig_count = zeros(M.ng,1);

for run_no = 1:P.permute_count
  fprintf('Run %d of %d\n', run_no, P.permute_count);
  tmp = randperm(M.np);
  P.patient_subset = tmp(1:P.permute_size);
  M = find_significant_genes(M,P);
  sig_genes = find(M.Q < P.Q_threshold);
  sig_count(sig_genes) = sig_count(sig_genes) + 1;
end

P_sig = sig_count / P.permute_count;
