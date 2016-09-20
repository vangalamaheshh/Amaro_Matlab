function M = take_patient_subset(M,pidx)

if islogical(pidx)
  pidx = find(pidx);
end
if ischar(pidx)
  pidx = {pidx};
end
if iscell(pidx)
  pidx = listmap(pidx,M.patient.name);
end
if isnumeric(pidx)
  if any(pidx<1 | pidx>slength(M.patient)), error('don''t understand pidx'); end
else
  error('don''t understand pidx');
end

pmap = listmap((1:slength(M.patient)),pidx);

M.patient = reorder_struct(M.patient,pidx);
M.np = length(pidx);

M.cov.sample = reorder_struct(M.cov.sample,pidx);
M.cov.ns = length(pidx);
M.cov.orig_cov = M.cov.orig_cov(pidx,:);
M.cov.fcov = M.cov.fcov(:,pidx);
M.cov.fcov_coding = M.cov.fcov_coding(:,pidx);
M.cov.gene_totcov = M.cov.gene_totcov(:,pidx);
M.cov.gene_cov = M.cov.gene_cov(:,pidx,:);
M.cov.gene_sil_cov = M.cov.gene_sil_cov(:,pidx,:);
M.cov.gene_non_cov = M.cov.gene_non_cov(:,pidx,:);
M.cov.gene_flank_cov = M.cov.gene_flank_cov(:,pidx,:);

M.mut.pat_idx = nansub(pmap,M.mut.pat_idx);

M.use(isnan(M.mut.pat_idx(M.use))) = [];
M.use_silent(isnan(M.mut.pat_idx(M.use_silent))) = [];
M.use_nonsilent(isnan(M.mut.pat_idx(M.use_nonsilent))) = [];
M.use_flank(isnan(M.mut.pat_idx(M.use_flank))) = [];

M.N_cov = M.N_cov(:,:,pidx);
M.N_sil_cov = M.N_sil_cov(:,:,pidx);
M.N_non_cov = M.N_sil_cov(:,:,pidx);
M.N_flank_cov = M.N_flank_cov(:,:,pidx);
M.n_flank = M.n_flank(:,:,pidx);
M.n_silent = M.n_silent(:,:,pidx);
M.n_nonsilent = M.n_nonsilent(:,:,pidx);
M.n_nonsilent_ignoring_null_categ = M.n_nonsilent_ignoring_null_categ(:,:,pidx);



