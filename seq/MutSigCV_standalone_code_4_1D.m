function G = MutSigCV_standalone_code_4_1D(M,gta,P,outdir)

% with 1D projection

% RESTORE THE USE OF BAGELS

flaghash = sparse([]);

% omitted gta?
if nargin==2 && (isempty(gta)||isstruct(gta)), P=gta; clear gta; end
if nargin==3 && (isempty(gta)||isstruct(gta)), outdir=P; P=gta; clear gta; end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'use_bagels_from',[]);
P = impose_default_value(P,'numbins_scheme',2);
P = impose_default_value(P,'use_flanking_mutations',false);
P = impose_default_value(P,'scale_flanking_mutations',true);

if exist('outdir','var')
  ede(outdir);
  progress_f = fopen([outdir '/progress.txt'],'wt');
  notes_f = fopen([outdir '/notes.txt'],'wt');
end

ng = M.ng;
np = M.np;
ncat = M.TOT-1;
ncols = ncat+1;

notes_fprintf('ng %d   np %d   ncat %d\n',ng,np,ncat);

if ~exist('gta','var'), gta = 1:ng; end
if ischar(gta), gta={gta}; end
if iscellstr(gta), gta = listmap(gta,M.gene.name); gta(isnan(gta))=[]; end

G=[];
G.gene = M.gene.name;
G.Nnon = round(sum(M.N_non_cov(:,end,:),3));
G.Nsil = round(sum(M.N_sil_cov(:,end,:),3));
G.Nflk = round(sum(M.N_flank_cov(:,end,:),3));
G.nnon = sum(M.n_nonsilent(:,end,:),3);
G.nsil = sum(M.n_silent(:,end,:),3);
G.nflk = sum(M.n_flank(:,end,:),3);

%%%%%%%%%%%%%%%
% TOTAL RATES %
%%%%%%%%%%%%%%%

globalrate_non = sum(G.nnon) / sum(G.Nnon);
globalrate_sil = sum(G.nsil) / sum(G.Nsil);
globalrate_flk = sum(G.nflk) / sum(G.Nflk);

notes_fprintf('Global rates (/Mb):  non %.2f   sil %.2f   flank %.2f\n',...
   globalrate_non*1e6,globalrate_sil*1e6,globalrate_flk*1e6);

%%%%%%%%%%%%%%
% BAGELS     %
%%%%%%%%%%%%%%

if ~isempty(P.use_bagels_from)
  fprintf('Loading existing bagels from %s\n',P.use_bagels_from);
  tmp = load(P.use_bagels_from);
  G.bagel = tmp.G.bagel;
  G.nnei = tmp.G.nnei;

else
  fprintf('Processing covariates...\n');
  nv = slength(M.V);
  V = nan(ng,nv);
  for vi=1:nv, V(:,vi) = M.V.val{vi}; end

  % convert covariate raw values to Z-scores
  Z = nan(ng,nv);
  for vi=1:nv
    missing = isnan(V(:,vi)) | isinf(V(:,vi));
    mn = mean(V(~missing,vi));
    sd = std(V(~missing,vi),0);  % second parameter=0 means confirm default behavior of normalize by (N-1) not (N)
    Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
  end
  
  fprintf('Finding bagels...  ');
keyboard

  min_neighbors = 0;
  max_neighbors = 50;
  qual_min = 0.05;
  
  G.nnei = zeros(ng,1); G.bagel = nan(ng,max_neighbors);
  for g=as_row(gta), if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
    
    % calculate distances from this gene
    df2 = bsxfun(@minus,Z,Z(g,:)).^2;
    dist2 = nansum(df2,2)./sum(~isnan(df2),2);
    [tmp,ord] = sort(dist2); ord = [g;ord(ord~=g)];
    
    % expand bagel outward until quality falls below qual_min
    nfit=0; Nfit=0;
    for ni=0:max_neighbors, gidx = ord(ni+1);
      
      Ngene = G.Nsil(gidx);
      ngene = G.nsil(gidx);
      if P.use_flanking_mutations
        ngene = ngene + G.nflk(gidx);
        if P.scale_flanking_mutations
          notes_fprintf_once('using silent (unscaled) + flanking (scaled) for bagels\n');
          Ngene = Ngene + round(G.Nflk(gidx)/(globalrate_non/globalrate_flk));
        else
          notes_fprintf_once('using silent (unscaled) + flanking (unscaled) for bagels\n');
          Ngene = Ngene + G.Nflk(gidx);
        end
      else
        notes_fprintf_once('using only silent (unscaled) for bagels\n');
      end
      
      if ni==0, ngene0=ngene; Ngene0=Ngene; end
      nfit=nfit+ngene; Nfit=Nfit+Ngene;
      
      % compare the gene being added to the central gene
      qual = 2*hyge2cdf(ngene,Ngene,ngene0,Ngene0);  % (two-sided)
      if qual>1, qual = 2-qual; end
      notes_fprintf_once('still using original two-tailed hyge2cdf as stopping criterion\n');
      
      % stopping criterion: stop if this gene would drop quality below qual_min
      if G.nnei(g)>=min_neighbors && qual<qual_min, break; end
      
      % update gene's bagel
      G.nnei(g) = ni;
      if ni>0
        G.bagel(g,ni)=gidx;
      end
      
    end % next neighborhood size
  end, fprintf('\n'); % next gene
  
end

% compute totals
fprintf('Totaling bagels:  ');
bagel_nsil = zeros(ng,ncat+1,np);
bagel_Nsil = zeros(ng,ncat+1,np);
bagel_nflk = zeros(ng,ncat+1,np);
bagel_Nflk = zeros(ng,ncat+1,np);
for g=as_row(gta), if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
  if G.nnei(g)>0
    gidx = G.bagel(g,1:G.nnei(g));
    bagel_nsil(g,:,:) = sum(M.n_silent(gidx,:,:),1);
    bagel_Nsil(g,:,:) = sum(M.N_sil_cov(gidx,:,:),1);
    bagel_nflk(g,:,:) = sum(M.n_flank(gidx,:,:),1);
    bagel_Nflk(g,:,:) = sum(M.N_flank_cov(gidx,:,:),1);
  end
end,fprintf('\n');

%%%%%%%%%%%%%%
% PROJECTION %
%%%%%%%%%%%%%%

notes_fprintf('Calculating p-value using 1D Projection method...\n');

%N_total = M.N_non_cov+M.N_sil_cov+M.N_flank_cov;
%n_total = M.n_nonsilent+M.n_silent+M.n_flank;
%notes_fprintf('Total = nonsilent+silent+flanking\n';

%n_total = M.n_nonsilent+M.n_silent;
%N_total = M.N_non_cov+M.N_sil_cov;
%notes_fprintf('Total:  n=nonsilent+silent   N=nonsilent+silent\n');

n_total = M.n_nonsilent+M.n_silent+bagel_nsil;
N_total = M.N_non_cov+M.N_sil_cov+bagel_Nsil;
if P.use_flanking_mutations
  ngene = ngene + (G.nflk+bagel_nflk);
  if P.scale_flanking_mutations
    N_total = N_total + (G.Nflk+G.bagel_Nflk)/(globalrate_non/globalrate_flk);
    notes_fprintf('Total:  n=nonsilent+silent+silent(bagel)+flk(bagel)  N=nonsilent+silent+silent(bagel)+flk(bagel)/scale\n');
  else
    N_total = N_total + (G.Nflk+G.bagel_Nflk);
    notes_fprintf('Total:  n=nonsilent+silent+silent(bagel)+flk(bagel)  N=nonsilent+silent+silent(bagel)+flk(bagel)\n');
  end
else
  notes_fprintf('Total:  n=nonsilent+silent+silent(bagel)   N=nonsilent+silent+silent(bagel)\n');
end

N_signal = M.N_non_cov;
n_signal = M.n_nonsilent;
notes_fprintf('Signal = nonsilent\n');

%N_signal = M.N_sil_cov;
%n_signal = M.n_silent;
%notes_fprintf('Signal = silent\n');

%N_signal = M.N_flank_cov;
%n_signal = M.n_flank;
%notes_fprintf('Signal = flanking\n');

% make sure values are rounded
N_total = round(N_total);
n_total = round(n_total);
N_signal = round(N_signal);
n_signal = round(n_signal);

% make sure we never have N>2*n
N_total(2*n_total>N_total)=2*n_total(2*n_total>N_total);
N_signal(2*n_signal>N_signal)=2*n_signal(2*n_signal>N_signal);
notes_fprintf('Making sure we never have N>2*n\n');

t1=0; t2=0; t3=0; t4=0;

G.pmin = nan(ng,1);
G.p = nan(ng,1);
G.pmax = nan(ng,1);
for g=as_row(gta) %, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end

  tt = tic;

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  P0 = nan(np,ncat);
  for c=1:ncat
    P0(:,c) = hygepdf(0,N_total(g,c,:),n_total(g,c,:),N_signal(g,c,:));
  end
  P0(isnan(P0))=1;
  P1 = 1-P0;

  % determine each patient's priority order of categories (according to P1)
  % left column of "priority" = least extreme category of mutation
  % right column of "priority" = most extreme category of mutation
  [tmp priority] = sort(P1,2,'descend');
  % sort the P arrays to put the columns in least->most priority order
  shft = (priority - repmat(1:ncat,np,1));
  map = reshape(1:(np*ncat),np,ncat);
  newmap = map + shft*np;
  P0 = P0(newmap);
  P1 = P1(newmap);

  t1=t1+toc(tt);
  tt = tic;

  % STEP 2
  % for each sample, compute probability that it would have been of each degree.
  % where d=0 (no mut) ..... ncat (most extreme mut)
  Pdeg = nan(np,ncat+1);
  for d=0:ncat
    % has to be clear in all categories > d
    Pdeg(:,d+1) = prod(P0(:,d+1:end),2);
    % and nonclear in category d (if d>0)
    if d>0, Pdeg(:,d+1) = Pdeg(:,d+1) .* P1(:,d); end
  end

  %% STEP 2a: calculate score for a sample being of each possible degree
  %% (uses new style, where score = -log10 probability of having a mutation in that category
  %% (zero score for no mutations)
  Sdeg = [zeros(np,1) -log10(P1)];

  % round all scores UP to nearest integer
  Sdeg = ceil(Sdeg);
  notes_fprintf_once('rounding all scores UP to nearest integer\n');

  % score cap
  Sdeg = min(Sdeg,1000);   % because inf will crash the C version!
  
  t2=t2+toc(tt);
  tt = tic;

  % STEP 3
  % determine actual degree of each sample, and score_obs for the gene

  degree = zeros(np,1);
  score_obs = 0;
  for p = 1:np
    for d = ncat:-1:1
      c = priority(p,d);
      if n_signal(g,c,p)>0, degree(p) = d; break; end
    end
    score_obs = score_obs + Sdeg(p,degree(p)+1);
  end

  % for zero score, don't bother doing convolutions
%  if score_obs<=0, G.p(g)=1; continue; end
%  notes_fprintf_once('skipping genes with zero mutations!\n');
  notes_fprintf_once('no longer skipping genes with zero mutations!\n');

  t3=t3+toc(tt);
  tt = tic;

  % STEP 4
  % compute P value for gene by convolutions
  % --> 95% of the total runtime is in this step (t4/ttot)

  if P.numbins_scheme==1
    numbins = max(10,score_obs*2);
    notes_fprintf_once('binsize fixed at 1; no minimum effect size.  numbins = max(10,score_obs*2)\n');
  elseif P.numbins_scheme==2
    numbins = ceil(score_obs+max(5,0.2*score_obs));
    notes_fprintf_once('binsize fixed at 1; no minimum effect size.  numbins = score_obs+max(5,0.2*score_obs)\n');
  else
    error('unknown P.numbins_scheme');
  end

  % allocate space for convolutions
  H = zeros(numbins,1);
  newH = zeros(numbins,ncols);

  [pmax pmin] = projection_1d_convolutions(Sdeg,Pdeg,score_obs,numbins,H,newH);
  % pmax is the p-value generated by the standard procedure we've always used.
  % pmin is the (better) p-value that goes one discrete step further
  % to solve the problem of discrete statistics, we want to take a randomly chosen position between pmin and pmax.
  G.pmax(g) = pmax;        % for sorting genelist
  G.p(g) = pmin+rand*(pmax-pmin);  % for q-q plot
  G.pmin(g) = pmin;        % for future reference

  t4=t4+toc(tt);
  tt = tic;
    
  ttot = t1+t2+t3+t4;
    
  progress_fprintf('%5d %12s   %8.6f %8.6f %8.6f   %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',...
          g,M.gene.name{g},G.pmin(g),G.p(g),G.pmax(g),t1,t2,t3,t4,ttot,ttot/g);

%keyboard

end, fprintf('\n');   % next gene

if exist('progress_f','var')
  fclose(progress_f);
  fclose(notes_f);
end

% FDR
G = sort_struct(G,'pmax');
G.q = calc_fdr_value(G.pmax);


% SAVE+RETURN RESULTS
if exist('outdir','var'), save([outdir '/results.mat'],'G'); end
G = rmfield(G,grep('bagel',fieldnames(G)));
if exist('outdir','var'), save_struct(G,[outdir '/sig_genes.txt']); end
pr(G,1:min(30,length(gta)))


  function progress_fprintf(varargin)
    fprintf(varargin{:});
    if exist('progress_f','var')
      fmt = regexprep(varargin{1},'(\S)\s+(\S)','$1\t$2');  % replace spaces with tabs
      fprintf(progress_f,fmt,varargin{2:end});
    end
  end
 
  function notes_fprintf(varargin)
    fprintf(varargin{:});
    if exist('notes_f','var'), fprintf(notes_f,varargin{:}); end
  end

  function notes_fprintf_once(varargin)
    hashval = sum(double(varargin{1}));
    if length(flaghash)>=hashval && flaghash(hashval)==1, return; end
    notes_fprintf(varargin{:});
    flaghash(hashval)=1;
  end

end
