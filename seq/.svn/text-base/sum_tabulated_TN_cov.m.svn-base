function c = sum_tabulated_TN_cov(pat,fsuf)
for i=1:slength(pat), fprintf('%d/%d ', i, slength(pat));
  tmp = load_matrix(['/xchip/tcga_scratch/lawrence/' pat.dir{i} '/' fsuf]);
  if i==1, c=tmp; else c=c+tmp; end
end,fprintf('\n');

