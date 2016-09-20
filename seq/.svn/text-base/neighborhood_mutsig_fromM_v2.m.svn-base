function M = neighborhood_mutsig_fromM_v2(M,P)

if ~exist('P','var'), P=[]; end

fprintf('NOTE: collapsing patients and categories for neighborhood_mutsig.\n');

if isfield(M.gene,'bagel')
  fprintf('Clearing existing bagels\n');
  M.gene = rmfield(M.gene,'bagel');
end

G=[]; G.rank=(1:slength(M.gene))'; G.name = M.gene.name;

np = slength(M.patient);
ng = slength(M.gene);

if isfield(M,'V')
  if length(M.V.val{1})~=ng, error('M.V and M.gene don''t agree in length'); end
else
  fprintf('Loading covariates data\n');
  tmp = load('/xchip/cga1/lawrence/mut/analysis/20110909_pancan/data.with_indels.v4.MVG.v2.mat','M','V');
  gidx = listmap(M.gene.name,tmp.M.gene.name);
  M.V = tmp.V;
  M.V = rmfield(M.V,'missing');
  for i=1:slength(M.V), M.V.val{i} = nansub(M.V.val{i},gidx); end
  clear tmp
end
V = M.V;

G.Nnon = sum(M.N_non_cov(:,end,:),3);
G.Nsil = sum(M.N_sil_cov(:,end,:),3);
G.Nflank = sum(M.N_flank_cov(:,end,:),3);

G.nnon = sum(M.n_nonsilent(:,end,:),3);
G.nsil = sum(M.n_silent(:,end,:),3);
G.nflank = sum(M.n_flank(:,end,:),3);


P = impose_default_value(P,'pval_calc_method','hyge2');
P = impose_default_value(P,'max_neighbors',30);
P = impose_default_value(P,'qual_min',0.05);
P = impose_default_value(P,'theta_halt',0.1);
P = impose_default_value(P,'incl_silent',true);
P = impose_default_value(P,'incl_flank',false);
P = impose_default_value(P,'incl_nonsilent',false);
P = impose_default_value(P,'theta_max',inf);
P = impose_default_value(P,'qual_metric','hyge2');
P = impose_default_value(P,'scale_for_negative_selection', false);
P = impose_default_value(P,'Z_zero_for_missing_data', false);
P = impose_default_value(P,'display_mode','mention');
P = impose_default_value(P,'p_mention_threshold', 0.01);

G = neighborhood_mutsig(G,V,P);

M.gene = mapinto(M.gene,G,'name');

