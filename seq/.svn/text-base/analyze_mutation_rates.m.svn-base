function M = analyze_mutation_rates(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'include_silent_in_sort',true);
P=impose_default_value(P,'exclude_these_genes_from_silent_calculation',{});
P=impose_default_value(P,'exclude_top_n_genes',40);
P=impose_default_value(P,'manual_correction_TSP_WashU',false);

fprintf('Analyzing mutation rates...\n');

genes_to_include = setdiff(1:M.ng,listmap(P.exclude_these_genes_from_silent_calculation,M.gene.name));
n_sil = sum(sum(M.n_silent(genes_to_include,M.TOT,:)));
N_tot = sum(sum(M.N_cov(genes_to_include,M.TOT,:)));
[mu_sil mu_sil_ci] = binofit(n_sil,N_tot);
sn_sil = sqrt(n_sil);
sN_tot = sqrt(N_tot);
mu_sil_sd = mu_sil * sqrt((sn_sil/n_sil)^2+(sN_tot/N_tot)^2);

fprintf('\nSilent mutation rate per total bases:\n\n');
fprintf('    n = %d\n', n_sil);
fprintf('    N = %d\n', N_tot);
fprintf('   mu = %5.2e +/- %1.2e  (95%%ci = %5.2e -- %5.2e)\n', ...
   mu_sil, mu_sil_sd, mu_sil_ci(1), mu_sil_ci(2));

% get relative rate calculations

C = compute_mutation_rates(M,P);

fprintf('\nRelative rates from silent mutations:\n\n');

n = C.Tnsil;
N = C.TNsil;
mu = n ./ N;
sn = n .^ 0.5;
sN = N .^ 0.5;
mu_sd = mu .* ((sn./n).^2+(sN./N).^2).^0.5;
mu_ci = 1.96 * mu_sd;
tot = repmat(mu(end),1,C.nc);
mu_rel = mu ./ tot;
mu_sd_rel = mu_sd ./ tot;
mu_ci_rel = 1.96 * mu_sd_rel;

fprintf('   class     n         N      sn    sN         mu +/- sd  (95%%ci)           mu_rel +/- sd  (95%%ci)\n');

for i=1:C.nc
  if i<C.nc, classname = M.mutclass{i};
  else classname = 'Total';
  end
  fprintf(['%8s   %3d %9d    %4.1f %5.0f   %6.2e +/- %6.2e (%6.2e) ' ...
           '%6.2f +/- %4.3f (%6.3f)\n'], ...
     classname, n(i), N(i), sn(i), sN(i), mu(i), mu_sd(i), mu_ci(i), mu_rel(i), mu_sd_rel(i), mu_ci_rel(i));
end

fprintf('\nRelative rates from non-top-%d genes plus silent mutations:\n\n', P.exclude_top_n_genes);

C = compute_mutation_rates(M,P);

n = C.npool;
N = C.Npool;
mu = n ./ N;
sn = n .^ 0.5;
sN = N .^ 0.5;
mu_sd = mu .* ((sn./n).^2+(sN./N).^2).^0.5;
mu_ci = 1.96 * mu_sd;
tot = repmat(mu(end),1,C.nc);
mu_rel = mu ./ tot;
mu_sd_rel = mu_sd ./ tot;
mu_ci_rel = 1.96 * mu_sd_rel;

fprintf('   class     n         N      sn    sN         mu +/- sd  (95%%ci)           mu_rel +/- sd  (95%%ci)\n');

for i=1:C.nc
  if i<C.nc, classname = M.mutclass{i};
  else classname = 'Total';
  end
  fprintf(['%8s   %3d %9d    %4.1f %5.0f   %6.2e +/- %6.2e (%6.2e) ' ...
           '%6.2f +/- %4.3f (%6.3f)\n'], ...
     classname, n(i), N(i), sn(i), sN(i), mu(i), mu_sd(i), mu_ci(i), mu_rel(i), mu_sd_rel(i), mu_ci_rel(i));
end

% print final calculation

fprintf('\nWeights applied to tallies of silent and nonsilent mutations from whole geneset:\n\n');

x = size(M.breakdown.frac,2);
tmp = round(sum(M.breakdown.frac .* repmat(M.N_terr(:,M.TOT),1,x)))';
c = [tmp(1:x/3) tmp(1+x/3:2*x/3)+tmp(1+2*x/3:end)];
cw = c .* repmat(mu_rel(1:end-1)',1,2);
cw_sd = c .* repmat(mu_sd_rel(1:end-1)',1,2);
ct = sum(cw);
ct_sd = sum(cw_sd.^2).^0.5;
ratio = ct(1)/ct(2);
ratio_sd = ratio * sqrt(sum((ct_sd./ct).^2));
mu_nsil = mu_sil / ratio;
mu_nsil_sd = mu_nsil * sqrt( (ratio_sd/ratio)^2 + (mu_sil_sd/mu_sil)^2);
mu_nsil_ci = 1.96 * mu_nsil_sd;
mu_nsil_high = mu_nsil + mu_nsil_ci;
mu_nsil_low = mu_nsil - mu_nsil_ci;

fprintf('   class    silent  nonsilent      weights (rel.rates)          silent               nonsilent\n');
for i=1:C.nc-1
  classname = M.mutclass{i};
  fprintf('%8s   %7d  %7d    x    %5.2f +/- %4.3f   =   %7.0f +/- %6.0f    %7.0f +/- %6.0f\n', ...
     classname, c(i,1),c(i,2),mu_rel(i),mu_sd_rel(i),cw(i,1),cw_sd(i,1),cw(i,2),cw_sd(i,2));
end

fprintf(['\n                                          ' ...
  '        TOTALS: %7.0f +/- %6.0f    %7.0f +/- %6.0f\n'], ct(1),ct_sd(1),ct(2),ct_sd(2));

fprintf(['\n                                          ' ...
  '    S/NS RATIO: %7.3f +/- %6.3f\n'], ratio, ratio_sd);

fprintf(['\n                                          ' ...
  '   silent rate: %7.2e +/- %6.2e\n'], mu_sil, mu_sil_sd);

fprintf(['\n                                          ' ...
  'nonsilent rate: %7.2e +/- %6.2e\n'], ...
  mu_nsil, mu_nsil_sd);

fprintf(['\n                                          ' ...
  '       (95%%ci = %4.2e -- %4.2e)\n\n'], ...
  mu_nsil_low, mu_nsil_high);


% save results of calculation into the data structure

% TOTAL RATE

M.mutrate.tot.low = mu_nsil_low;
M.mutrate.tot.hat = mu_nsil;
M.mutrate.tot.high = mu_nsil_high;
M.mutrate.tot.stdev = mu_nsil_sd;


% RELATIVE RATES

%   for non-indel classes

rel = C.RRpool(1:C.nc-1,C.HAT)';

%   for indel classes

n_indel = sum(sum(M.n_indel(C.nontop,M.TOT-M.NUM_INDEL_CLASSES:M.TOT-1,:),3),1);
N_indel = repmat(sum(sum(M.N_cov(C.nontop,M.TOT,:),3),1),1,M.NUM_INDEL_CLASSES);

indel_tot = n_indel ./ N_indel;
indel_rel = indel_tot / mu_nsil;

rel = [rel indel_rel];

%   for total

rel = [rel 1.0];

M.mutrate.rel = rel;
