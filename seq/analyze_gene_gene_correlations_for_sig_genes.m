function R = analyze_gene_gene_correlations_for_sig_genes(M,P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'mutsig_significance_threshold',0.1);
P = impose_default_value(P,'max_genes_to_consider',100);

if ~isfield(M.gene,'qval')
  error('Requires significance analysis to have been run');
end

sig_genes = find(M.gene.qval<=P.mutsig_significance_threshold);

if length(sig_genes)<2
  fprintf('Need at least 2 significant genes in order to perform correlation analysis!\n');
  R = [];
  return;
end

if length(sig_genes)>P.max_genes_to_consider
  fprintf('Truncating list of significant genes after the first %d genes\n',P.max_genes_to_consider);
  [tmp ord] = sort(M.gene.qval);
  sig_genes = ord(1:P.max_genes_to_consider);
end

P.genes_to_consider = M.gene.name(sig_genes);
R = analyze_gene_gene_correlations_4(M,P);

