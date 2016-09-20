function G = MutSigCV_standalone_code_2(M,gta,P,outdir)
% changed API to use M
% M = new_load_mutdata(maf,P);

% omitted gta?
if nargin==2 && (isempty(gta)||isstruct(gta)), P=gta; clear gta; end
if nargin==3 && (isempty(gta)||isstruct(gta)), outdir=P; P=gta; clear gta; end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'use_s_directly_instead_of_bagels',false);
P = impose_default_value(P,'use_sf_directly_instead_of_bagels',false);
P = impose_default_value(P,'use_ns_directly_instead_of_bagels',false);
P = impose_default_value(P,'use_nsf_directly_instead_of_bagels',false);
P = impose_default_value(P,'use_precalculated_bagels','');
P = impose_default_value(P,'max_neighbors',50);
P = impose_default_value(P,'min_neighbors',0);
P = impose_default_value(P,'scale_silent_and_flanking_rates',false);
P = impose_default_value(P,'scale_silent_and_flanking_rates_by_adjusting_N',false);
P = impose_default_value(P,'additional_scale_silent',1);
P = impose_default_value(P,'additional_scale_flanking',1);
P = impose_default_value(P,'qual_min',0.05);
P = impose_default_value(P,'set_Xgsc_equal_to_Ngsc',false);
P = impose_default_value(P,'blurfactor',1);
P = impose_default_value(P,'projection_bin_score_cap',inf);
P = impose_default_value(P,'projection_probability_method','hyge2');
P = impose_default_value(P,'null_score_boost',3);
P = impose_default_value(P,'min_effect_size',1.25); if P.min_effect_size<1, error('P.min_effect_size set incorrectly'); end
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

globalrate_non = sum(G.nnon) / sum(G.Nnon);
globalrate_sil = sum(G.nsil) / sum(G.Nsil);
globalrate_flk = sum(G.nflk) / sum(G.Nflk);
globalrate_sf = (sum(G.nflk+G.nsil)) / (sum(G.Nflk+G.Nsil));

% What background model to use?  (direct counts or bagels)

if P.use_s_directly_instead_of_bagels+P.use_sf_directly_instead_of_bagels+...
      P.use_nsf_directly_instead_of_bagels+P.use_ns_directly_instead_of_bagels>1
  error('please choose only one');

elseif P.use_s_directly_instead_of_bagels
  fprintf('Using silent directly instead of bagels\n');
  if P.scale_silent_and_flanking_rates_by_adjusting_N
    fprintf('                (with scaling of N)\n');
    x_gsc = M.n_silent;
    X_gsc = M.N_sil_cov/(globalrate_non/globalrate_sil*P.additional_scale_silent);
  elseif P.scale_silent_and_flanking_rates
    fprintf('                (with scaling of counts)\n');
    x_gsc = M.n_silent*(globalrate_non/globalrate_sil*P.additional_scale_silent);
    X_gsc = M.N_sil_cov;
  else
    x_gsc = M.n_silent;
    X_gsc = M.N_sil_cov;
  end

elseif P.use_sf_directly_instead_of_bagels
  fprintf('Using silent+flanking directly instead of bagels\n');
  if P.scale_silent_and_flanking_rates_by_adjusting_N
    fprintf('                (with scaling of N)\n');
    x_gsc = M.n_silent + M.n_flank;
    X_gsc = M.N_sil_cov/(globalrate_non/globalrate_sil*P.additional_scale_silent) + ...
            M.N_flank_cov/(globalrate_non/globalrate_flk*P.additional_scale_flanking);
  elseif P.scale_silent_and_flanking_rates
    fprintf('                (with scaling of counts)\n');
    x_gsc = M.n_silent*(globalrate_non/globalrate_sil*P.additional_scale_silent) + ...
            M.n_flank*(globalrate_non/globalrate_flk*P.additional_scale_flanking);
    X_gsc = M.N_sil_cov + M.N_flank_cov;
  else
    x_gsc = M.n_silent + M.n_flank;
    X_gsc = M.N_sil_cov + M.N_flank_cov;
  end

else  % use bagels
  if ~isempty(P.use_precalculated_bagels)
    fprintf('Loading precalculated bagels from %s\n',P.use_precalculated_bagels)
    tmp = load(P.use_precalculated_bagels,'G');
    G.nnei = tmp.G.nnei;
    G.nfit = tmp.G.nfit;
    G.Nfit = tmp.G.Nfit;
  else
    fprintf('Processing covariates...\n');

    nv = slength(M.V);
    V = nan(ng,nv);
    for vi=1:nv
      if isfield(M.V,'val2')
        G.(M.V.name{vi}) = M.V.val2{vi};   % original (untransformed) value for display in table
      else
        G.(M.V.name{vi}) = M.V.val{vi};
      end
      V(:,vi) = M.V.val{vi};
    end

    % cosmetic changes to covariates for improved display by pr()
    if isfield(G,'smoothed_CCLE_expression')
      G = rename_field(G,'smoothed_CCLE_expression','expr');
      G.expr = round(G.expr);
    end
    if isfield(G,'rt_extra1_max')
      G = rename_field(G,'rt_extra1_max','RT');
    end
    if isfield(G,'hiC')
      G = make_numeric(G,'hiC');
      if all(abs(G.hiC)<1), G.hiC = round(G.hiC*1000); end
    end

    % convert covariate raw values to Z-scores
    Z = nan(ng,nv);
    for vi=1:nv
      missing = isnan(V(:,vi)) | isinf(V(:,vi));
      mn = mean(V(~missing,vi));
      sd = std(V(~missing,vi),0);  % second parameter=0 means confirm default behavior of normalize by (N-1) not (N)
      Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
    end

    fprintf('Finding bagels...  ');

    max_neighbors = P.max_neighbors;
    min_neighbors = P.min_neighbors;
    qual_min = P.qual_min;
    
    G.nnei = zeros(ng,1); G.nfit = zeros(ng,1); G.Nfit = zeros(ng,1);
    
    for g=as_row(gta), if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
      
      % calculate distances from this gene
      df2 = bsxfun(@minus,Z,Z(g,:)).^2;
      dist2 = nansum(df2,2)./sum(~isnan(df2),2);
      [tmp,ord] = sort(dist2); ord = [g;ord(ord~=g)];
      
      % expand bagel outward until quality falls below qual_min
      nfit=0; Nfit=0;
      for ni=0:max_neighbors, gidx = ord(ni+1);
        
        Ngene = G.Nsil(gidx) + G.Nflk(gidx);
        if P.scale_silent_and_flanking_rates
          ngene = G.nsil(gidx)*(globalrate_non/globalrate_sil) + ...
                  G.nflk(gidx)*(globalrate_non/globalrate_flk);
        else
          ngene = G.nsil(gidx) + G.nflk(gidx);
        end
        if ni==0, ngene0=ngene; Ngene0=Ngene; end
        nfit=nfit+ngene; Nfit=Nfit+Ngene;
        
        % compare the gene being added to the central gene
        qual = 2*hyge2cdf(ngene,Ngene,ngene0,Ngene0);  % (two-sided)
        if qual>1, qual = 2-qual; end
        
        % stopping criterion: stop if this gene would drop quality below qual_min
        if G.nnei(g)>=min_neighbors && qual<qual_min, break; end
        
        % update gene's statistics
        G.nnei(g) = ni; G.nfit(g) = nfit; G.Nfit(g) = Nfit;
        
      end % next neighborhood size
    end, fprintf('\n'); % next gene
  end

  fprintf('Expanding to (x,X)_gsc...\n');

  n_gsc = M.n_nonsilent + M.n_silent + M.n_flank;
  N_gsc = M.N_non_cov + M.N_sil_cov + M.N_flank_cov;
  n_sc = sum(n_gsc,1);
  N_sc = sum(N_gsc,1);
  n_c = sum(n_sc,3);
  N_c = sum(N_sc,3);
  mu_c = n_c./N_c;
  n_tot = n_c(end);
  N_tot = N_c(end);
  mu_tot = n_tot/N_tot;
  f_c = mu_c/mu_tot;
  f_Nc = N_c/N_tot;
  n_s = n_sc(:,end,:);
  N_s = N_sc(:,end,:);
  mu_s = n_s./N_s;
  f_s = mu_s/mu_tot;
  f_Ns = N_s/mean(N_s);
  x_gsc = repmat(G.nfit,[1 ncat+1 np]); X_gsc = repmat(G.Nfit,[1 ncat+1 np]);       % last column = total
  x_gsc = bsxfun(@times,x_gsc,f_c.*f_Nc); X_gsc = bsxfun(@times,X_gsc,f_Nc);
  x_gsc = bsxfun(@times,x_gsc,f_s.*f_Ns); X_gsc = bsxfun(@times,X_gsc,f_Ns);
  
  if P.set_Xgsc_equal_to_Ngsc
    % BEFORE:
    % round(X_gsc(1:10,:,1))
    %      250188      499982     2447458     4894867     2848216     5696334    16636515    16636515
    %       97428      194703      953087     1906155     1109150     2218261     6478577     6478577
    %      133800      267390     1308897     2617768     1523222     3046392     8897185     8897185
    %      583213     1165508     5705270    11410425     6639476    13278725    38781381    38781381
    %      450197      899685     4404039     8807990     5125175    10250176    29936306    29936306
    %       14479       28934      141637      283270      164829      329652      962770      962770
    %       87838      175538      859275     1718532      999976     1999918     5840891     5840891
    %      122109      244026     1194531     2389037     1390128     2780209     8119781     8119781
    %      169035      337805     1653586     3307139     1924352     3848638    11240197    11240197
    %      127701      255201     1249232     2498438     1453786     2907523     8491609     8491609
    
    fprintf('      (setting Xgsc=Ngsc)\n');
    mu_gsc = x_gsc ./ X_gsc;
    X_gsc = M.N_sil_cov + M.N_flank_cov;
    x_gsc = mu_gsc .* X_gsc;
    
    % AFTER:
    % round(X_gsc(1:10,:,1))
    %          19          32         105         127          36          52         369         369
    %           8          17         112         140         175         229         682         682
    %          23          40          92         119         132         207         612         612
    %          15          22         301         347         307         433        1423        1423
    %          22          41         359         422         349         503        1697        1697
    %          15          22          76          80          11          15         219         219
    %           4           6          74          70          51          62         267         267
    %          10          20         190         260         145         230         854         854
    %          25          38         185         233         113         168         762         762
    %           4           7          57          64         109         128         370         370
    % ratios are ~20000-40000 !
    % why so high? would expect them to be ~np  (5300) i guess it's because they're multiplied by both np and nnei
    
    %---> result of using this method: all genes get p=1
  end
end

%%%%% add F columns
G.Fnon = (G.nnon./G.Nnon) / globalrate_non;
G.Fsil = (G.nsil./G.Nsil) / globalrate_sil;
G.Fflk = (G.nflk./G.Nflk) / globalrate_flk;
if isfield(G,'nfit') && isfield(G,'Nfit')
  if P.scale_silent_and_flanking_rates
    G.Ffit = (G.nfit./G.Nfit) / globalrate_non;
  else
    G.Ffit = (G.nfit./G.Nfit) / globalrate_sf;
  end
end

%%%%%%%%%%%%%%
% PROJECTION %
%%%%%%%%%%%%%%

fprintf('Calculating p-value using 2D Projection method...\n');

blurfactor = P.blurfactor;
null_score_boost = P.null_score_boost;
min_effect_size = P.min_effect_size;
convolution_numbins = P.convolution_numbins;

if ncat==1 && null_score_boost~=0
  fprintf('Only one category: no null boost\n');
  null_score_boost = 0;
end

x_gsc = x_gsc / blurfactor;
X_gsc = X_gsc / blurfactor;

G.p = nan(ng,1);

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

  N = reshape(M.N_non_cov(g,1:ncat,:),ncat,np)';
  n = reshape(M.n_nonsilent(g,1:ncat,:),ncat,np)';
  N(2*n>N)=2*n(2*n>N);  % make sure we don't have N>2*n

  x = reshape(x_gsc(g,1:ncat,:),ncat,np)';
  X = reshape(X_gsc(g,1:ncat,:),ncat,np)';

  if strcmp(P.projection_probability_method,'hyge2')
    P0 = hyge2pdf(0,N,x,X);
    P1 = hyge2pdf(1,N,x,X);

  elseif strcmp(P.projection_probability_method,'hyge')
    P0 = nan(np,ncat);
    P1 = nan(np,ncat);
    for p=1:np
      for c=1:ncat
        nx = round(n(p,c)+x(p,c));
        NX = round(N(p,c)+X(p,c));
        Npc = round(N(p,c));
        if nx==0
          P0(p,c) = 1;
          P1(p,c) = 0;
        else
          P0(p,c) = hygepdf(0,NX,nx,Npc);
          P1(p,c) = hygepdf(1,NX,nx,Npc);
          if isnan(P0(p,c))||isnan(P1(p,c))
            fprintf('what?');
            P0(p,c) = 1;
            P1(p,c) = 0;
          end
        end
      end
    end        

  else
    error('unknown P.projection_probability_method');
  end

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

  % null score boost
  % NOTE: has a problem, only adds null boost for d1=null, d2=0
  priority2 = [zeros(np,1) priority];
  Sdeg(priority2==ncat) = Sdeg(priority2==ncat) + null_score_boost;

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

%keyboard

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








