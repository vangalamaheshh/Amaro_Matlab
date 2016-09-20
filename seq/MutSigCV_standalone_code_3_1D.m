function G = MutSigCV_standalone_code_3_1D(M,gta,P,outdir)

% with 1D projection

% omitted gta?
if nargin==2 && (isempty(gta)||isstruct(gta)), P=gta; clear gta; end
if nargin==3 && (isempty(gta)||isstruct(gta)), outdir=P; P=gta; clear gta; end

if ~exist('P','var'), P=[]; end

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
% PROJECTION %
%%%%%%%%%%%%%%

notes_fprintf('Calculating p-value using 1D Projection method...\n');

G.pmin = nan(ng,1);
G.pmax = nan(ng,1);
G.p = nan(ng,1);

%N_total = M.N_non_cov+M.N_sil_cov+M.N_flank_cov;
%n_total = M.n_nonsilent+M.n_silent+M.n_flank;
%notes_fprintf('Total = nonsilent+silent+flanking\n';

n_total = M.n_nonsilent+M.n_silent;
N_total = M.N_non_cov+M.N_sil_cov;
notes_fprintf('Total:  n=nonsilent+silent   N=nonsilent+silent\n');

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

G.p = nan(ng,1);
firstgene=true;
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
  shft = (priority - repmat([1:ncat],np,1));
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
  if firstgene
    notes_fprintf('  rounding all scores UP to nearest integer\n');
  end

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
  if score_obs<=0, G.p(g)=1; continue; end
  if firstgene
    notes_fprintf('skipping genes with zero mutations!\n');
%    notes_fprintf('no longer skipping genes with zero mutations!\n');
  end

  t3=t3+toc(tt);
  tt = tic;

  % STEP 4
  % compute P value for gene by convolutions
  % --> 95% of the total runtime is in this step (t4/ttot)

  if firstgene
    notes_fprintf('binsize fixed at 1; no minimum effect size.  numbins = max(10,score_obs*2)\n');
  end
  numbins = max(10,score_obs*2);

  % allocate space for convolutions
  H = zeros(numbins,1);
  newH = zeros(numbins,ncols);

  [pmax pmin] = projection_1d_convolutions(Sdeg,Pdeg,score_obs,numbins,H,newH);
  % p1 is the p-value generated by the standard procedure we've always used.
  % p2 is the (better) p-value that omits the density of the "=obs" bin.
  % to solve the problem of discrete statistics, we want to take a randomly chosen position between p1 and p2.
  randfrac = rand;
  p = pmin+randfrac*(pmax-pmin);
  
  G.pmax(g) = pmax;
  G.pmin(g) = pmin;
  G.p(g) = p;

  t4=t4+toc(tt);
  tt = tic;
    
  ttot = t1+t2+t3+t4;
    
  progress_fprintf('%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n',...
          g,M.gene.name{g},G.p(g),G.pmax(g),G.pmin(g),t1,t2,t3,t4,ttot,ttot/g);

%keyboard

  firstgene=false;
end, fprintf('\n');   % next gene

if exist('progress_f','var')
  fclose(progress_f);
  fclose(notes_f);
end

% FDR
G.q = calc_fdr_value(G.p);
G = sort_struct(G,'p');

if exist('outdir','var')
  save([outdir '/results.mat'],'G');
  save_struct(G,[outdir '/sig_genes.txt']);
end

pr(G,1:min(30,length(gta)))



  function progress_fprintf(varargin)
    fprintf(varargin{:});
    if exist('progress_f','var') fprintf(progress_f,varargin{:}); end
  end
 
  function notes_fprintf(varargin)
    fprintf(varargin{:});
    if exist('notes_f','var') fprintf(notes_f,varargin{:}); end
  end

end
