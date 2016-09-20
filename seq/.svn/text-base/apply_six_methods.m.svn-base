function apply_six_methods(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'filename_stem','sig');

P.sig_calculation_method = 'concatenation';
rate = {'low';'hat';'high'};

P.use_1D_concatenation = true;
for r=1:length(rate)
  P.mutation_rate_to_use = rate{r};
  M = find_significant_genes(M,P);
  P.filename = [P.filename_stem '_simple_' rate{r} '.txt'];
  output_sig_table(M,P);
end

P.use_1D_concatenation = false;
for r=1:length(rate)
  P.mutation_rate_to_use = rate{r};
  M = find_significant_genes(M,P);
  P.filename = [P.filename_stem '_categ_' rate{r} '.txt'];
  output_sig_table(M,P);
end

end
