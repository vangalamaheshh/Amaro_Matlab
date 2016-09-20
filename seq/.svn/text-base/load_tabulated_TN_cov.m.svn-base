function c = load_tabulated_TN_cov(pat,fsuf)
for i=1:slength(pat), fprintf('%d/%d ', i, slength(pat));
  c{i,1} = load_matrix(['/xchip/tcga_scratch/lawrence/' pat.dir{i} '/' fsuf]);
end,fprintf('\n');

