function G = MutSigCV_standalone_code(mutation_file,coverage_file,covariate_file,output_file)

if nargin~=4, error('usage: MutSigCV_standalone_code(mutation_file,coverage_file,covariate_file,output_file)'); end

fprintf('Loading data...\n');

M = load_struct(mutation_file);
M = make_numeric(M,'categ');
M = make_boolean(M,{'is_coding','is_silent'});
G = load_struct_specify_string_cols(covariate_file,1);
C = load_struct_specify_string_cols(coverage_file,1:2);

fprintf('Building n and N tables...\n');

f = fieldnames(C); patient_names = f(4:end);
f = fieldnames(G); cvnames = f(2:end);

ng = slength(G);
np = length(patient_names);
ncat = max(M.categ);
nv = length(cvnames);

% make sure C is sorted by the same gene order as in G
C.gene_idx = listmap(C.gene,G.gene);
C = sort_struct(C,'gene_idx');

M.gene_idx = listmap(M.Hugo_Symbol,G.gene);
M.patient = regexprep(M.Tumor_Sample_Barcode,'-Tumor$','');
M.patient = regexprep(M.patient,'-','_');
M.patient_idx = listmap(M.patient,patient_names);

midx = find(M.is_coding & M.is_silent);
n_silent = hist3d(M.gene_idx(midx),M.categ(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
midx = find(M.is_coding & ~M.is_silent);
n_nonsilent = hist3d(M.gene_idx(midx),M.categ(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
midx = find(~M.is_coding);
n_flank = hist3d(M.gene_idx(midx),M.categ(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);

N_silent = nan(ng,ncat,np);
N_nonsilent = nan(ng,ncat,np);
N_flank = nan(ng,ncat,np);

for ci=1:ncat
  silent_idx = find(strcmp(C.zone,'silent') & C.categ==ci);
  nonsilent_idx = find(strcmp(C.zone,'nonsilent') & C.categ==ci);
  flank_idx = find(strcmp(C.zone,'flank') & C.categ==ci);
  for pi=1:np
    N_silent(:,ci,pi) = C.(patient_names{pi})(silent_idx);
    N_nonsilent(:,ci,pi) = C.(patient_names{pi})(nonsilent_idx);
    N_flank(:,ci,pi) = C.(patient_names{pi})(flank_idx);
  end
end

% add total columns
n_silent(:,end+1,:) = sum(n_silent,2);
n_nonsilent(:,end+1,:) = sum(n_nonsilent,2);
n_flank(:,end+1,:) = sum(n_flank,2);
N_silent(:,end+1,:) = N_silent(:,end,:);          % copy total coverage from null coverage
N_nonsilent(:,end+1,:) = N_nonsilent(:,end,:);
N_flank(:,end+1,:) = N_flank(:,end,:);

% total across patients, save in G
G.N_nonsilent = sum(N_nonsilent(:,end,:),3);
G.N_silent = sum(N_silent(:,end,:),3);
G.N_flank = sum(N_flank(:,end,:),3);
G.n_nonsilent = sum(n_nonsilent(:,end,:),3);
G.n_silent = sum(n_silent(:,end,:),3);
G.n_flank = sum(n_flank(:,end,:),3);

fprintf('Processing covariates...\n');

V = nan(ng,nv);
for vi=1:nv, V(:,vi) = G.(cvnames{vi}); end

% convert covariate raw values to Z-scores
Z = nan(ng,nv);
for vi=1:nv
  missing = isnan(V(:,vi)) | isinf(V(:,vi));
  mn = mean(V(~missing,vi));
  sd = std(V(~missing,vi),0);  % second parameter=0 confirms the default of normalizing by (N-1) not (N)
  Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
end

fprintf('Finding bagels...  ');

max_neighbors = 50;
qual_min = 0.05;

G.nnei = nan(ng,1); G.x = nan(ng,1); G.X = nan(ng,1);

for gi=1:ng, if ~mod(gi,1000), fprintf('%d/%d ',gi,ng); end

  % calculate distances from this gene
  df2 = bsxfun(@minus,Z,Z(gi,:)).^2;
  dist2 = nansum(df2,2)./sum(~isnan(df2),2);
  [tmp,ord] = sort(dist2); ord = [gi;ord(ord~=gi)];

  % expand bagel outward until quality falls below qual_min
  nfit=0; Nfit=0;
  for ni=0:max_neighbors, gidx = ord(ni+1);

    ngene = G.n_silent(gidx) + G.n_flank(gidx);
    Ngene = G.N_silent(gidx) + G.N_flank(gidx);
    if ni==0, ngene0=ngene; Ngene0=Ngene; end
    nfit=nfit+ngene; Nfit=Nfit+Ngene;

    % compare the gene being added to the central gene
    hc = hyge2cdf(ngene,Ngene,ngene0,Ngene0);
    qual_left = min(hc, 1-hc);
    qual = 2*qual_left;

    % stopping criterion: stop if this gene would drop quality below qual_min
    if ni>0 && qual<qual_min, break; end

    % update gene's statistics
    G.nnei(gi) = ni; G.x(gi) = nfit; G.X(gi) = Nfit;

  end % next neighborhood size
end, fprintf('\n'); % next gene

fprintf('Expanding to (x,X)_gcp...\n');

n_gcp = n_nonsilent + n_silent + n_flank;
N_gcp = N_nonsilent + N_silent + N_flank;

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

x_gcp = repmat(G.x,[1 ncat+1 np]); X_gcp = repmat(G.X,[1 ncat+1 np]);       % last column = total
x_gcp = bsxfun(@times,x_gcp,f_c.*f_Nc); X_gcp = bsxfun(@times,X_gcp,f_Nc);
x_gcp = bsxfun(@times,x_gcp,f_p.*f_Np); X_gcp = bsxfun(@times,X_gcp,f_Np);

fprintf('Calculating p-value using 2D Projection method...  ');

null_score_boost = 3;
min_effect_size = 1.25;
convolution_numbins = 1000;

G.p = nan(ng,1);

for g=1:ng, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  N = reshape(N_nonsilent(g,1:ncat,:),ncat,np)';
  n = reshape(n_nonsilent(g,1:ncat,:),ncat,np)';
  N(n>N)=n(n>N);  % make sure we don't have N>n

  x = reshape(x_gcp(g,1:ncat,:),ncat,np)';
  X = reshape(X_gcp(g,1:ncat,:),ncat,np)';
  P0 = hyge2pdf(0,N,x,X);
  P1 = hyge2pdf(1,N,x,X);

  % determine each patient's priority order of categories (according to P1)
  % left column of "priority" = least extreme category of mutation
  % right column of "priority" = most extreme category of mutation
  [tmp priority] = sort(P1,2,'descend');
  % sort the P arrays to put the columns in least->most priority order
  shft = (priority - repmat([1:ncat],np,1));
  map = reshape(1:(np*ncat),np,ncat);
  newmap = map + shft*np;
  P0 = P0(newmap);
  P1 = P1(newmap);
  P2 = 1-(P0+P1);  % note, P2 means "P(2 or more)"
  P2(P2<0) = 0;

  % STEP 2
  % for each sample, compute probability that it would have been of each (2-dimensional) degree.
  % degree=(d1,d2), where d=0 (no mut) ..... ncat (most extreme mut)
  % d1 is the MOST extreme mutation (or no mutation)
  % d2 is the SECOND MOST extreme mutation (or no mutation)
  % d1 can be 0-ncat; d2 can be 0-d1

  Pdeg = zeros(np,ncat+1,ncat+1);
  for d1=0:ncat, for d2=0:d1
    % has to have 0 in any/all categories > d1
    p = prod(P0(:,d1+1:end),2);
    if (d1>0)  % and (if d1>0)
      if (d1==d2)
        % if d1==d2, has to have 2+ in category d1
        p = p .* P2(:,d1);
      else
        % else:   has to have exactly 1 in category d1
        %         has to be clear in any/all categories (d2+1) to (d1-1)
        %         and (if d2>0) have (1 or 2+) in category d2
        p = p .* P1(:,d1);
        p = p .* prod(P0(:,d2+1:d1-1),2);
        if (d2>0)
          p = p .* (P1(:,d2)+P2(:,d2));
        end
      end
    end
    Pdeg(:,d1+1,d2+1) = p;
    end,end

  %% STEP 2a: calculate score for a sample being of each possible degree
  %% (uses new style, where score = -log10 probability of having a mutation in that category
  %% (zero score for no mutations)
  Sdeg = zeros(np,ncat+1,ncat+1);
  for d1=1:ncat, for d2=0:d1
    if d1==d2
      p = P2(:,d1);
    else
      if d2>0
        p = P1(:,d1).*P1(:,d2);
      else
        p = P1(:,d1);
      end
    end
    Sdeg(:,d1+1,d2+1) = -log10(p);
  end,end

  % null score boost
  priority2 = [zeros(np,1) priority];
  Sdeg(priority2==ncat) = Sdeg(priority2==ncat) + null_score_boost;

 % STEP 3
  % determine actual (two-dimensional) degree and score for each sample
  % sum scores to get score_obs for gene

  degree = zeros(np,2);
  score_obs = 0;
  for p = 1:np
    i = 1;
    for d = ncat:-1:1
      c = priority(p,d);
      if i==1
        if n(p,c)>=2
          degree(p,:) = [d d];
          i = 3;
        elseif n(p,c)==1
          degree(p,i) = d;
          i=i+1;
        end
      elseif i==2
        if n(p,c)>=1
          degree(p,i) = d;
          i=i+1;
        end
      else % i>2: done
        break
      end
    end
    score_sample = Sdeg(p,degree(p,1)+1,degree(p,2)+1);
    score_obs = score_obs + score_sample;
  end

  % impose minimum effect size by decreasing score_obs
  score_obs = score_obs / min_effect_size;

  % for zero score, don't bother doing convolutions
  if score_obs<=0, G.p(g)=1; continue; end

  % STEP 4
  % compute P value for gene by convolutions

  numbins = convolution_numbins;
  binsize = score_obs / numbins;
  H = zeros(numbins,1);
  H(1) = 1;  % initial condition: all probability is in first bin

  % sequential convolution
  offset = min(numbins, round(Sdeg/binsize));
  ncols = (ncat+1)*(ncat+2)/2;
  newH = zeros(numbins,ncols);
  for p=1:np
    newH(:) = 0;
    col=1;
    for d1=0:ncat, for d2=0:d1
      o = offset(p,d1+1,d2+1);
      newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
      col=col+1;
    end,end
    H = sum(newH,2);
  end

  % save p-value
  G.p(g) = max(0,1-sum(H));

end, fprintf('\n');   % next gene

% FDR
G.q = calc_fdr_value(G.p);

G = sort_struct(G,'p');
save_struct(G,output_file);

fprintf('Done.\n');








