function G = MutSigCV_standalone_code_2(M,gta,P,outdir)

% omitted gta?
if nargin==2 && (isempty(gta)||isstruct(gta)), P=gta; clear gta; end
if nargin==3 && (isempty(gta)||isstruct(gta)), outdir=P; P=gta; clear gta; end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'projection_bin_score_cap',inf);
P = impose_default_value(P,'convolution_numbins',1000);

if exist('outdir','var')
  ede(outdir);
  out_f = fopen([outdir '/progress.txt'],'wt');
end

ng = M.ng;
np = M.np;
ncat = M.TOT-1;

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

globalrate_non = sum(G.nnon) / sum(G.Nnon)
globalrate_sil = sum(G.nsil) / sum(G.Nsil)
globalrate_flk = sum(G.nflk) / sum(G.Nflk)

%%%%%%%%%%%%%%
% PROJECTION %
%%%%%%%%%%%%%%

fprintf('Calculating p-value using 2D Projection method...\n');

convolution_numbins = P.convolution_numbins;

G.p = nan(ng,1);

%N_total = round(M.N_non_cov+M.N_sil_cov+M.N_flank_cov);
%n_total = round(M.n_nonsilent+M.n_silent+M.n_flank);
%N_signal = round(M.N_non_cov);
%n_signal = round(M.n_nonsilent);

N_total = round(M.N_non_cov+M.N_sil_cov);
n_total = round(M.n_nonsilent+M.n_silent);
fprintf('Total = nonsilent+silent\n');

%N_signal = round(M.N_non_cov);
%n_signal = round(M.n_nonsilent);

N_signal = round(M.N_sil_cov);
n_signal = round(M.n_silent);
fprintf('Signal = silent\n');

%N_signal = round(M.N_flank_cov);
%n_signal = round(M.n_flank);

% make sure we never have N>2*n
N_total(2*n_total>N_total)=2*n_total(2*n_total>N_total);
N_signal(2*n_signal>N_signal)=2*n_signal(2*n_signal>N_signal);

t1=0; t2=0; t3=0; t4=0;

% allocate space for convolutions
H = zeros(convolution_numbins,1);
ncols = (ncat+1)*(ncat+2)/2;
newH = zeros(convolution_numbins,ncols);

for g=as_row(gta) %, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end

  tt = tic;

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  P0 = nan(np,ncat);
  P1 = nan(np,ncat);
  for c=1:ncat
    P0(:,c) = hygepdf(0,N_total(g,c,:),n_total(g,c,:),N_signal(g,c,:));
    P1(:,c) = hygepdf(1,N_total(g,c,:),n_total(g,c,:),N_signal(g,c,:));
  end
  P0(isnan(P0))=1;
  P1(isnan(P1))=0;

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

  t1=t1+toc(tt);
  tt = tic;

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

  % score cap
  Sdeg = min(Sdeg,P.projection_bin_score_cap);

  t2=t2+toc(tt);
  tt = tic;

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
        if n_signal(g,c,p)>=2
          degree(p,:) = [d d];
          i = 3;
        elseif n_signal(g,c,p)==1
          degree(p,i) = d;
          i=i+1;
        end
      elseif i==2
        if n_signal(g,c,p)>=1
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

  % for zero score, don't bother doing convolutions
  if score_obs<=0, G.p(g)=1; continue; end

  t3=t3+toc(tt);
  tt = tic;

  % STEP 4
  % compute P value for gene by convolutions
  % --> 95% of the total runtime is in this step (t4/ttot)

  Sdeg(isinf(Sdeg)) = 1000;  % to stop the C version from crashing!
  G.p(g) = projection_convolutions_fast(Sdeg,Pdeg,score_obs,convolution_numbins,H,newH);
  if G.p(g)<0, G.p(g)=0; end

  t4=t4+toc(tt);
  tt = tic;

  ttot = t1+t2+t3+t4;

  fff ={'%d\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n',...
          g,M.gene.name{g},G.p(g),t1,t2,t3,t4,ttot,ttot/g};
  fprintf(fff{:}); if exist('out_f','var'), fprintf(out_f,fff{:}); end

end, fprintf('\n');   % next gene

if exist('out_f','var')
  fclose(out_f);
end

% FDR
G.q = calc_fdr_value(G.p);
G = sort_struct(G,'p');

if exist('outdir','var')
  save([outdir '/results.mat'],'G');
end

pr(G,1:min(30,length(gta)))








