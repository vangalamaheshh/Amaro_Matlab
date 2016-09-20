function hmutnames = identify_hypermutated_samples(M,P)

% default parameters

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'exclude_these_genes',{});
P=impose_default_value(P,'number_of_IQRs',2);

fprintf('Identifying hypermutated samples\n');

hist_genes = 1:M.ng;
exclude_list = listmap(P.exclude_these_genes, M.gene.name);
hist_genes = setdiff(hist_genes, exclude_list);

N_breakdown = repmat(M.N_cov(:,M.TOT,:),1,24) .* repmat(M.breakdown.frac, [1 1 M.np]);
N_si = round(squeeze(sum(sum(N_breakdown(hist_genes,1:8,:),2))));
N_ns = round(squeeze(sum(sum(N_breakdown(hist_genes,9:24,:),2))));

n_ns = squeeze(sum(M.n_nonsilent(hist_genes,M.TOT,:),1));
n_si = squeeze(sum(M.n_silent(hist_genes,M.TOT,:),1));

rtot = (n_ns + n_si) ./ (N_ns + N_si);
hmut = outliers(rtot, struct('method','iqr','n',P.number_of_IQRs));
hmutnames = M.patient.name(hmut);
