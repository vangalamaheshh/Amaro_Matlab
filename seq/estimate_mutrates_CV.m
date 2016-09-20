function M = estimate_mutrates_CV(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'key','value');
P = impose_default_value(P,'patients_to_analyze',1:slength(M.patient));
P = impose_default_value(P,'use_corrected_expansion_method',true);

if isfield(M,'mutrate_analysis')
  fprintf('Overwriting previous mutrate_analysis\n');
  M = rmfield(M,'mutrate_analysis');
end

if ~isfield(M,'n_flank')
  M.n_flank = M.n_nonsilent;
  M.n_flank(:) = 0;
end
if ~isfield(M,'N_flank_cov')
  M.N_flank_cov = M.N_non_cov;
  M.N_flank_cov(:) = 0;
end

pidx = P.patients_to_analyze;

M.gene.x1 = repmat({'---------------------'},slength(M.gene),1);

M.gene.Nflank = sum(M.N_flank_cov(:,end,pidx),3);
M.gene.Nsil = sum(M.N_sil_cov(:,end,pidx),3);
M.gene.Nnon = sum(M.N_non_cov(:,end,pidx),3);
M.gene.nflank = sum(M.n_flank(:,end,pidx),3);
M.gene.nsil = sum(M.n_silent(:,end,pidx),3);
M.gene.nnon = sum(M.n_nonsilent(:,end,pidx),3);

M.gene.x2 = repmat({'---------------------'},slength(M.gene),1);


% calculate nfit and Nfit using the bagels
fprintf('Calculating nfit and Nfit\n');
nfit = nan(M.ng,1); Nfit = nan(M.ng,1);
for i=1:M.ng, if ~mod(i,1000), fprintf('%d/%d ',i,M.ng); end
  b = [i M.gene.bagel(i,:)]; b(isnan(b)) = [];
  nfit(i) = fullsum(M.n_silent(b,end,:)) + fullsum(M.n_flank(b,end,:));
  Nfit(i) = fullsum(M.N_sil_cov(b,end,:)) + fullsum(M.N_flank_cov(b,end,:));
end, fprintf('\n');
M.gene.Nfit = Nfit;
M.gene.nfit = nfit;
n_tot = fullsum(M.n_silent(:,end,pidx)) + fullsum(M.n_flank(:,end,pidx));
N_tot = fullsum(M.N_sil_cov(:,end,pidx)) + fullsum(M.N_flank_cov(:,end,pidx));
mu_tot = sum(n_tot)/sum(N_tot);
M.gene.fMLE = (M.gene.nfit./M.gene.Nfit)/mu_tot; % (for easy reference)

%% check for problem cases:
%% (1) nnon>0, Nfit==0       --> should be taken care of when MutSig removes poorly covered genes
idx = find(M.gene.nnon>0 & M.gene.Nfit==0);
if ~isempty(idx), fprintf('PROBLEM: some genes have nnon>0 and Nfit==0\n'); keyboard; error('please fix'); end
idx = find(M.gene.nnon==0 & M.gene.Nfit==0);
if ~isempty(idx), fprintf('WARNING: some genes have Nfit==0 (but also nnon==0)\n'); keyboard; end

% calculate per-category and per-sample rates
fprintf('Calculating x_gcp,X_gcp\n');
nc = M.TOT-1;
ng = size(M.N_sil_cov,1);
np = size(M.N_sil_cov,3);

if ~P.use_corrected_expansion_method

  % original expansion method
  fprintf('MutSigCV: WARNING: Using original expansion method, which has semicolon bug!\n');

  % (has SEMICOLON BUG!)                                     _
  n = M.n_nonsilent + M.n_silent + M.n_flank; N = M.N_non_cov; + M.N_sil_cov + M.N_flank_cov;

  n = sum(n,1); n = n(:,1:nc,pidx); N = sum(N,1); N = N(:,1:nc,pidx);
  mu_tot = fullsum(n)./fullsum(N); mu_c = sum(n,3)./sum(N,3); mu_s = sum(n,2)./sum(N,2);
  f_c = mu_c./mu_tot; f_s = mu_s./mu_tot;
  f_s(f_s==0) = min(f_s(f_s>0));  % for samples with zero mutations, replace with the lowest nonzero rate
  X_c = sum(N,3); X_s = sum(N,2); fX_c = X_c/mean(X_c); fX_s = X_s/mean(X_s);

  % calculate final priors
  x_g = M.gene.nfit; X_g = M.gene.Nfit;
  x = bsxfun(@times,bsxfun(@times,x_g,f_c.*fX_c),f_s.*fX_s); X = bsxfun(@times,bsxfun(@times,X_g,fX_c),fX_s);
  
else

  % corrected expansion method
  % ("METHOD6")
  fprintf('MutSigCV: Using corrected expansion method.\n');

  n_gcp = M.n_nonsilent + M.n_silent + M.n_flank;
  N_gcp = M.N_non_cov + M.N_sil_cov + M.N_flank_cov;
  
  n_cp = sum(n_gcp,1);
  N_cp = sum(N_gcp,1);

  n_c = sum(n_cp,3);
  N_c = sum(N_cp,3);
  mu_c = n_c./N_c;

  n_tot = n_c(end);
  N_tot = N_c(end);
  mu_tot = n_tot/N_tot;
  f_c = mu_c/mu_tot;
  f_Nc = N_c/N_tot;

  n_p = n_cp(:,end,:);
  N_p = N_cp(:,end,:);
  mu_p = n_p./N_p;
  f_p = mu_p/mu_tot;
  f_Np = N_p/mean(N_p);

  x = repmat(M.gene.nfit,[1 nc+1 np]); X = repmat(M.gene.Nfit,[1 nc+1 np]);       % last column = total
  x = bsxfun(@times,x,f_c.*f_Nc); X = bsxfun(@times,X,f_Nc);
  x = bsxfun(@times,x,f_p.*f_Np); X = bsxfun(@times,X,f_Np);

  % remove total column
  x(:,nc+1,:) = [];
  X(:,nc+1,:) = [];

end

% fix weird problem with x>X
idx = find(x(:)>X(:));
if ~isempty(idx)
  fprintf('PROBLEM: %d/%d entries have x>X\n',length(idx),numel(x));
  fprintf('Will handle this by increasing X to 1000x in these cases\n');
  X(x(:)>X(:)) = 1000 * x(x(:)>X(:));
%  keyboard; error('please fix'); end
end

M.mutrate_analysis = [];
M.mutrate_analysis.x = x;
M.mutrate_analysis.X = X;



