function M = get_mutsig_pvals(M,gta)

if ~exist('gta','var')
  gta = M.gene.name;
  gi = (1:slength(M.gene))';
else
  if isnumeric(gta)
    gta = M.g.name(1:gta);
  elseif ischar(gta)
    gta = {gta};
  end
  gi = listmap(gta,M.gene.name);
  gi(isnan(gi)) = [];
end

P=[];
P.use_prior = true;
P.patient_names = M.sample.name;
P.gene_names = gta;
ng = length(gta);
np = slength(M.sample);
ncat = M.ncat-1;   % (exclude Total)

% data about category-specific mutation rate
M.categ.ntot = histc(M.mut.categ(~isnan(M.mut.gene_idx)),1:slength(M.categ));
M.categ.ntot(end) = sum(M.categ.ntot);
M.categ.Ntot = squeeze(sum(sum(M.gene_sil_cov + M.gene_non_cov + M.gene_flank_cov,1),2));
M.categ.mu = (M.categ.ntot./M.categ.Ntot);
M.categ.murel = M.categ.mu/M.categ.mu(end);
M.categ.murel(isnan(M.categ.murel)) = 1;

% data about sample-specific mutation rate
M.sample.ntot = histc(M.mut.pat_idx(~isnan(M.mut.gene_idx)),1:slength(M.sample));
M.sample.Ntot = sum(M.gene_sil_cov(:,:,end) + M.gene_non_cov(:,:,end) + M.gene_flank_cov(:,:,end),1)';
M.sample.mu = (M.sample.ntot./M.sample.Ntot);
mu_tot = (sum(M.sample.ntot)/sum(M.sample.Ntot));
M.sample.murel = M.sample.mu/mu_tot;
M.sample.murel(isnan(M.sample.murel)) = 1;

% product of sample- and category-specific marginals
murel_catpat = bsxfun(@times,M.categ.murel(1:ncat)',M.sample.murel);

% data about gene-specific rate

%% for each gene, choose background model that gave the *highest* conservative estimate
nbkgd = [M.gene.nbagelflank_c(gi,end) M.gene.nbagelsil_c(gi,end) M.gene.nflank_c(gi,end) ...
    M.gene.nsil_c(gi,end) M.gene.nsf_c(gi,end) M.gene.nsphere_c(gi,end)];
Nbkgd = [M.gene.Nbagelflank_c(gi,end) M.gene.Nbagelsil_c(gi,end) M.gene.Nflank_c(gi,end) ...
    M.gene.Nsil_c(gi,end) M.gene.Nsf_c(gi,end) M.gene.Nsphere_c(gi,end)];
fbkgd = [M.gene.fbagelflank_c(gi,end) M.gene.fbagelsil_c(gi,end) M.gene.fflank_c(gi,end) ...
    M.gene.fsil_c(gi,end) M.gene.fsf_c(gi,end) M.gene.fsphere_c(gi,end)];
[tmp model] = max(fbkgd,[],2);

x=nan(ng,1); X=nan(ng,1);
for i=1:ng
  x(i) = nbkgd(i,model(i));
  X(i) = Nbkgd(i,model(i));
end
x(X==0) = 0;

% actual n and N nonsilent counts
N=nan(ng,ncat,np); n=nan(ng,ncat,np);
new_gene_idx = listmap(M.mut.gene_idx,gi);
for c=1:ncat
  idx = find(M.mut.categ==c & M.mut.is_coding & ~M.mut.is_silent);
  n(:,c,:) = hist2d_fast(new_gene_idx(idx),M.mut.pat_idx(idx),1,ng,1,np);
  N(:,c,:) = M.gene_non_cov(gi,:,c);
end
n(N==0) = 0;

% add prior from what the per-gene rates are like
if P.use_prior
%  [M.prior.a,M.prior.b]=mle_beta(M.gene.x,M.gene.X);
%  fprintf('Priors:  a = %f   b = %f\n',M.prior.a,M.prior.b);
%  M.gene.x_sphere = M.gene.x_sphere + (M.prior.a-1);
%  M.gene.X_sphere = M.gene.X_sphere + (M.prior.b+M.prior.a-2);
end


%% scale by mu_cs
X = repmat(X,[1 ncat np]);
x = bsxfun(@times,x,shiftdim(murel_catpat',-1));

%%%%%% call MutSig

M.gene.mutsig_p = nan(slength(M.gene),1);
M.gene.mutsig_q = nan(slength(M.gene),1);

N = round(N);
n = round(n);

keff = 1.5;

p = calculate_significance(N,n,{X,keff*x},P);

if ng<50
  pr(P.gene_names,p);
end


if nargout>0
  M.gene.mutsig_p(gi) = p;
  M.gene.mutsig_q(gi) = calc_fdr_value(p);
  flds = {'mutsig_p','mutsig_q'};
  M.g = mapinto(M.g,M.gene,'name',flds);
  M.g = orderfields_last(M.g,flds);
  M.gene = orderfields_last(M.gene,flds);
end


