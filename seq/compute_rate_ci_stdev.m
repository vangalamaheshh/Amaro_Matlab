function F = compute_rate_ci_stdev(F)
[r c] = binofit_2d(F.n,F.N); F.rate = 1e6*r; F.cilow = 1e6*c(:,:,1); F.cihigh = 1e6*c(:,:,2); F.stdev = (F.cihigh-F.rate)/1.96;
F.effect = bsxfun(@rdivide,F.rate,F.rate(1,:));
