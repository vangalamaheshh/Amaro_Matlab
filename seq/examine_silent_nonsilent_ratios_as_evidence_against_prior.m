function M = examine_silent_nonsilent_ratios_as_evidence_against_prior(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'impute_full_coverage',~isfield(M,'N_cov'));

fprintf('examine_silent_nonsilent_ratios_as_evidence_against_prior... ');

M.mutrate.per_gene_BMR_correction = ones(M.ng,1);

global_non_bmr = M.mutrate.tot.hat;


x=[];
x.exp = (-8:0.01:-2)';
x.bmr = 10.^x.exp;
%prior_strength = 1e3;

%for g=grep('RYR2',M.gene.name,1)
for g=1:M.ng, if ~mod(g,1000), fprintf('%d/%d ',g,M.ng);end
  if ~P.impute_full_coverage
    Ncov = fullsum(M.N_cov(g,M.TOT,:));
  else
    Ncov fullsum(M.N_terr).*M.np;
  end  
  nsil = fullsum(M.n_silent(g,M.TOT,:));
  s_ns_ratio = fullsum(M.N_sil_cov(g,M.TOT,:))./fullsum(M.N_non_cov(g,M.TOT,:));
  global_sil_bmr = global_non_bmr*s_ns_ratio;
x.prior = normpdf(x.bmr,global_sil_bmr,1e-3);
%  x.prior = betapdf(x.bmr,prior_strength*global_sil_bmr,prior_strength);
  x.like = binopdf(nsil,Ncov,x.bmr);
  x.tot = x.like.*x.prior;
  idx = find(x.bmr>=global_sil_bmr,1);
  x.tot(1:idx-1) = 0;
  [tmp idx] = max(x.tot);
  M.mutrate.per_gene_BMR_correction(g) = x.bmr(idx)/global_sil_bmr;
% M.mutrate.per_gene_BMR_correction(g)
end, fprintf('\n');

