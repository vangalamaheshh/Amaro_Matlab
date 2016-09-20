function N = load_and_sum_covlines(samples,covfile)
% N = load_and_sum_covlines(samples)
%
% Mike Lawrence 2010-01-27

for i=1:length(samples)
  fn = ['/xchip/tcga_scratch/lawrence/' samples{i} '/' covfile];
  tmp = load_matrix(fn);
  tmp = tmp(6:end);  % get rid of patient gene chr start end
  if i==1, N=tmp; else N=N+tmp; end
end
N = N';
