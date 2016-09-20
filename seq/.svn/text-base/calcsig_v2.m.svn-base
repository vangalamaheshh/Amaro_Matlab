function M = calcsig_v2(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'use_bagel',true);
P = impose_default_value(P,'use_prior',true);
P = impose_default_value(P,'scale_magnitudes',true);
P = impose_default_value(P,'use_category_specific_BMRs',true);
P = impose_default_value(P,'use_sample_specific_BMRs',true);
P = impose_default_value(P,'take_max_of_gene_bagel_sphere',false);
P = impose_default_value(P,'take_max_of_gene_bagel_sphere_and_scale_down_X',false);
P = impose_default_value(P,'collapse_categories',false);
P = impose_default_value(P,'collapse_samples',false);

ng = slength(M.gene);
ncat = M.TOT-1;
np = slength(M.patient);

% actual nonsilent mutations and territory
M.gene.N_work = M.N_non_cov(:,1:ncat,:);
M.gene.n_work = M.n_nonsilent(:,1:ncat,:);

if P.collapse_categories
  M.gene.N_work = M.gene.N_work(:,end,:);
  M.gene.n_work = sum(M.gene.n_work,2);
  ncat = 1;
end

if P.collapse_samples
  M.gene.N_work = sum(M.gene.N_work,3);
  M.gene.n_work = sum(M.gene.n_work,3);
  np = 1;
end

% data about gene-specific mutation rate
M.gene.x = M.gene.nsil + M.gene.nflank;
M.gene.X = M.gene.Nsil + M.gene.Nflank;
[M.gene.rate M.gene.ci_low M.gene.ci_high] = binofit_2d(M.gene.x,round(max(M.gene.x,M.gene.X)));
globalrate = sum(M.gene.x)/sum(M.gene.X);
M.gene.f = M.gene.ci_low/globalrate;

M.gene.x_bagel = max(0,M.gene.nfit - M.gene.x);
M.gene.X_bagel = max(0,M.gene.Nfit - M.gene.X);

if P.take_max_of_gene_bagel_sphere

  q = {'gene'   'bagel'        'sphere'}';
  x = [M.gene.x M.gene.x_bagel M.gene.nfit];
  X = [M.gene.X M.gene.X_bagel M.gene.Nfit];
  r = x./X;
  r(X==0) = nan;

  [tmp idx] = max(r,[],2);

  M.gene.bkgd_type = nansub(q,idx);
  M.gene.x_work = nan(ng,1);
  M.gene.X_work = nan(ng,1);
  for i=1:ng
    M.gene.x_work(i) = x(i,idx(i));
    M.gene.X_work(i) = X(i,idx(i));
  end

  if P.take_max_of_gene_bagel_sphere_and_scale_down_X
    [tmp idx2] = min(X,[],2);
    fac = nan(ng,1);
    for i=1:ng
      fac(i) = X(idx2(i))/X(idx(i));
    end
    M.gene.x_work = M.gene.x_work .* fac;
    M.gene.X_work = M.gene.X_work .* fac;
  end

elseif P.use_bagel

  % compute dependence of "w" upon "f"
  [f_sort ord] = sort(M.gene.f);
  binsize = 200;   % genes per bin
  nbins = ceil(ng/binsize);
  div = unique(f_sort(1:binsize:end));
  bin=[];bin.min = div; bin.max = [div(2:end);max(M.gene.f)+1];
  nbins = slength(bin);
  ww = 0:0.01:1;
  l=nan(length(ww),1);
  M.gene.w_bagel = nan(ng,1);
  for i=1:nbins
    gidx = find(M.gene.f>=bin.min(i) & M.gene.f<bin.max(i));
    bin.gidx{i,1} = gidx;
    bin.ngenes(i,1) = length(gidx);
    for j=1:length(l)
      l(j)=bbinologlik(M.gene.x(gidx),M.gene.X(gidx),ww(j)*M.gene.x_bagel(gidx)+1,...
              ww(j)*(M.gene.X_bagel(gidx)-M.gene.x_bagel(gidx))+1);
    end
    [tmp idx] = max(l);
    bin.w(i,1) = ww(idx);
    M.gene.w_bagel(gidx) = ww(idx);
  end

  % combine gene + bagel into a single x+X:
  M.gene.x_work = M.gene.x + M.gene.w_bagel.*M.gene.x_bagel;
  M.gene.X_work = M.gene.X + M.gene.w_bagel.*M.gene.X_bagel;

else % no bagel

  M.gene.x_work = M.gene.x;
  M.gene.X_work = M.gene.X;
  
end

% add prior from what the per-gene rates are like
if P.use_prior
  [M.prior.a,M.prior.b]=mle_beta(M.gene.x,M.gene.X);
  fprintf('Priors:  a = %f   b = %f\n',M.prior.a,M.prior.b);
  M.gene.x_work = M.gene.x_work + (M.prior.a-1);
  M.gene.X_work = M.gene.X_work + (M.prior.b+M.prior.a-2);
end

M.gene.X_work = repmat(M.gene.X_work,[1 ncat np]);
M.gene.x_work = repmat(M.gene.x_work,[1 ncat np]);

% adjustments to x,X

if P.scale_magnitudes || P.use_category_specific_BMRs
  % data about category-specific mutation rate
  % (based on nonsilent)
  M.categ.ntot = sum(sum(M.n_nonsilent(:,1:ncat,:),3),1)';
  M.categ.Ntot = sum(sum(M.N_non_cov(:,1:ncat,:),3),1)';
  M.categ.mu = (M.categ.ntot./M.categ.Ntot);
  n_tot = sum(M.categ.ntot);
  N_tot = M.categ.Ntot(end);
  mu_tot = n_tot/N_tot;
  M.categ.Nrel = M.categ.Ntot/N_tot;
  M.categ.murel = M.categ.mu/mu_tot;

  if P.scale_magnitudes
    % scale magnitude to the per-category territories
    M.gene.X_work = bsxfun(@times,M.gene.X_work,M.categ.Nrel');
    M.gene.x_work = bsxfun(@times,M.gene.x_work,M.categ.Nrel');
  end

  % if requested, scale rate (numerator) to the per-category rates
  if P.use_category_specific_BMRs
    M.gene.x_work = bsxfun(@times,M.gene.x_work,M.categ.murel');
  end
end

if P.use_sample_specific_BMRs
  % data about sample-specific mutation rate
  % (based on nonsilent)
  M.patient.ntot = squeeze(sum(M.n_nonsilent(:,end,:),1));
  M.patient.Ntot = squeeze(sum(M.N_non_cov(:,end,:),1));
  M.patient.mu = (M.patient.ntot./M.patient.Ntot);
  mu_tot = (sum(M.patient.ntot)/sum(M.patient.Ntot));
  M.patient.murel = M.patient.mu/mu_tot;
  % if requested, also scale rate (numerator) to the per-sample rates
  M.gene.x_work = bsxfun(@times,M.gene.x_work,shiftdim(M.patient.murel,-2));
end

% call MutSig
P.gene_names = M.gene.name;
M.gene.ntot = sum(M.gene.n_work,3);
M.gene.Ntot = sum(M.gene.N_work(:,end,:),3);
M.gene = rmfield_if_exist(M.gene,{'p','q','rank'});
M.gene = orderfields_first(M.gene,{'name'});
if isfield(M.gene,'X_work') && isfield(M.gene,'x_work') && ~isfield(M.gene,'mu_work')
  M.gene.p = calculate_significance(M.gene.N_work,M.gene.n_work,{M.gene.X_work,M.gene.x_work},P);
elseif ~isfield(M.gene,'X_work') && ~isfield(M.gene,'x_work') && isfield(M.gene,'mu_work')
  M.gene.p = calculate_significance(M.gene.N_work,M.gene.n_work,M.gene.mu_work,P);
else
  error('M.gene needs either mu_work OR X_work, x_work');
end

M.gene.q = calc_fdr_value(M.gene.p);

if isfield(M.gene,'effect')
  [tmp ord] = sort_struct(M.gene,{'p','effect'},[1 -1]);
else
  [tmp ord] = sort_struct(M.gene,'p');
end

[tmp M.gene.rank] = sort(ord);


if isfield(M,'g')
  M.g = mapinto(M.g,M.gene,'name',{'p','q'},{'p2','q2'})
end
