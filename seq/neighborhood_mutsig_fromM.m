function M = neighborhood_mutsig_fromM(M,varargin)

fprintf('NOTE: collapsing patients and categories for neighborhood_mutsig.\n');

G=[]; G.rank=(1:slength(M.gene))'; G.name = M.gene.name;

np = slength(M.pat);
ng = slength(M.gene);

if isfield(M,'V')
  V = M.V;
  if length(V.val{1})~=ng, error('M.V and M.gene don''t agree in length'); end
else
  fprintf('Loading covariates data\n');
  tmp = load('/xchip/cga1/lawrence/mut/analysis/20110909_pancan/data.with_indels.v4.MVG.v2.mat','M','V');
  gidx = listmap(M.gene.name,tmp.M.gene.name);
  V = tmp.V;
  V = rmfield(V,'missing');
  for i=1:slength(V), V.val{i} = nansub(V.val{i},gidx); end
  clear tmp
end

G.Nnon = M.cov.gene_non_repcov(:,end)*np;
G.Nsil = M.cov.gene_sil_repcov(:,end)*np;
G.Nflank = M.cov.gene_flank_repcov(:,end)*np;

G.nnon = histc(M.mut.gene_idx(M.mut.is_coding & ~M.mut.is_silent),1:ng);
G.nsil = histc(M.mut.gene_idx(M.mut.is_coding & M.mut.is_silent),1:ng);
G.nflank = histc(M.mut.gene_idx(~M.mut.is_coding & M.mut.is_flank),1:ng);

G = neighborhood_mutsig(G,V,varargin{:});

M.gene = mapinto(M.gene,G,'name');

