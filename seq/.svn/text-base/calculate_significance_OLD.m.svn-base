function M = calculate_significance(M,P)

if nargout>1
  error('output format of calculate_significance has been changed!');
end

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'sig_calculation_method','concatenation');
P=impose_default_value(P,'use_1D_concatenation',false);
P=impose_default_value(P,'use_semi_exact_procedure',true);
P=impose_default_value(P,'semi_exact_procedure_bin_size',0.01);
P=impose_default_value(P,'mutation_rate_to_use','hat');
P=impose_default_value(P,'manual_mutation_rate',NaN);
P=impose_default_value(P,'zero_top_genes',false);
P=impose_default_value(P,'top_genes_to_zero',{});
P=impose_default_value(P,'patient_subset',1:M.np);
P=impose_default_value(P,'use_sample_specific_mutation_rates',false);
P=impose_default_value(P,'make_sample_specific_mutation_rates_equal_to_global',false);
P=impose_default_value(P,'convolution_maxmuts',1000);
P=impose_default_value(P,'convolution_maxmuts_limit',100000);
P=impose_default_value(P,'first_gene_to_calculate',1);
P=impose_default_value(P,'last_gene_to_calculate',M.ng);
P=impose_default_value(P,'pval_cutoff',1e-11);

multimode = isfield(M,'multi');
if multimode
  if strcmp('sig_calculation_method','projection'), error('projection not implemented for multimode'); end
  if length(unique(P.patient_subset))~=M.np, error('patient subsets not implemented for multimode'); end
  nsets = length(M.multi);
  check_multiM_agreement(M.multi);
end

fprintf('Computing significance of each gene\n');

% BUILD WORKING TABLES

if ~multimode
  nmutcats = M.TOT-1;
  N_work = M.N_cov(:,1:nmutcats,P.patient_subset);
  n_work = M.n_nonsilent(:,1:nmutcats,P.patient_subset);
else % multimode
  nmutcats = 0;
  for i=1:nsets
    nmutcats = nmutcats + (M.multi{i}.TOT-1);
  end
  N_work = zeros(M.ng,nmutcats,M.np);
  n_work = N_work;
  coffset = 0;
  poffset = 0;
  for i=1:nsets
    this_np = M.multi{i}.np;
    this_ncat = M.multi{i}.TOT-1;
    N_work(:,coffset+[1:this_ncat],poffset+[1:this_np]) = M.multi{i}.N_cov(:,1:this_ncat,:);
    n_work(:,coffset+[1:this_ncat],poffset+[1:this_np]) = M.multi{i}.n_nonsilent(:,1:this_ncat,:);
    coffset = coffset + this_ncat;
    poffset = poffset + this_np;
  end
end

if P.use_1D_concatenation
  Ntot_cov = sum(M.N_cov(:,:,P.patient_subset),3);
  ntot_nonsilent = sum(M.n_nonsilent(:,:,P.patient_subset),3);
end

% CHOOSE MUTATION RATES

demand_field(M,'mutrate_analysis');
if P.use_sample_specific_mutation_rates
  if multimode
    error('sample-specific rates not implemented with multimode');
  end
  if ~strcmpi(P.sig_calculation_method,'concatenation') || P.use_1D_concatenation==true
    error('Sample-specific rates implemented only with concatenation-with-categories method');
  end
  if P.make_sample_specific_mutation_rates_equal_to_global
    fprintf('WARNING: setting all sample-specific mutation rates equal to the global rate\n');
    rate = repmat(sum(M.mutrate_analysis.ss.n,1)./sum(M.mutrate_analysis.ss.N,1),M.np,1);
  else
    rate = M.mutrate_analysis.ss.rate;
  end
  % check to make sure we won't need to take log(0)
  tmp_tot = squeeze(sum(n_work,1))';
  [a b] = find(tmp_tot>0 & rate==0);
  if ~isempty(a), error('Zero mutation rate in bins with nonzero mutation counts'); end
  if isfield(M,'exprcorr')
    error('Expression-based BMR correction not yet implemented with sample-specific mutation rates');
  end 
else
  if strcmpi(P.sig_calculation_method,'Fisher')
    if multimode, error('Fisher method not implemented with multimode'); end
    demand_field(M.mutrate_analysis,{'n','N'});
    rate_n = M.mutrate_analysis.n;
    rate_N = M.mutrate_analysis.N;
  else
    if ~multimode
      demand_field(M.mutrate_analysis,{'tot','rel'});
      if strcmp(P.mutation_rate_to_use,'hat')
        rate_tot = M.mutrate_analysis.tot.hat;
      elseif strcmp(P.mutation_rate_to_use,'high')
        rate_tot = M.mutrate_analysis.tot.high;
      elseif strcmp(P.mutation_rate_to_use,'low')
        rate_tot = M.mutrate_analysis.tot.low;
      elseif strcmp(P.mutation_rate_to_use,'manual')
        if isnan(P.manual_mutation_rate), error('Must specify manually-set mutation rate'); end
        rate_tot = P.manual_mutation_rate;
      else
        error('Please specify mutation_rate_to_use = high/hat/low/manual');
      end
      rate = rate_tot * M.mutrate_analysis.rel(1:nmutcats);
    else % multimode
      if ~strcmp(P.mutation_rate_to_use,'hat'), error('non-hat rates not implemented with multimode'); end
      rate = cell(nsets,1);
      for i=1:nsets, rate{i} = M.multi{i}.mutrate_analysis.rate; end
      rate = cat(2,rate{:});
      rate_tot = M.mutrate_analysis.tot.hat;
    end
    % check to make sure we won't need to take log(0)
    tmp_tot = sum(sum(n_work,3),1);
    for c=1:nmutcats
      if tmp_tot(c)>0 && rate(c)==0
        fprintf('--> mutation category %s, with %d mutations, has rate of zero!', M.mutclass{c}, tmp(c));
        fprintf('Can correct this problem by setting mu_excludes_known_cancer_genes = false\n');
        fprintf('  or by manually setting indel rate.\n');
        error('Zero mutation rate');
      end
    end
  end
  % EXPRESSION-BASED BMR CORRECTIONS
  if isfield(M,'exprcorr')
    if strcmpi(P.sig_calculation_method,'Fisher')
      fprintf('WARNING: Ignoring expression-based BMR correction in Fisher method\n');
    else
      BMR_adjustment_table = M.exprcorr;
      if size(BMR_adjustment_table,1) ~= M.ng, error('M.exprcorr should have one row per gene'); end
      if size(BMR_adjustment_table,2) ~= length(rate), error('M.exprcorr should have one column per category'); end
    end
  else
    BMR_adjustment_table = ones(M.ng,length(rate));
  end
  % expand rate from (category) to (gene x category)
  rate = bsxfun(@times,rate,BMR_adjustment_table);
end

if isfield(M.mutrate_analysis,'per_gene_BMR_correction')
  fprintf('Applying final per-gene BMR correction\n');
  rate = bsxfun(@times,rate,M.mutrate_analysis.per_gene_BMR_correction);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% for the rest of the function, only n_work, N_work, and rate (or rate_tot [1D], or rate_n,rate_N [Fisher]) are used.
%
%

% select appropriate method and begin

Prob = nan(M.ng,1);
s_g = nan(M.ng,1);

if strcmpi(P.sig_calculation_method, 'projection')
  METHOD = 1;
elseif strcmp(P.sig_calculation_method, 'Fisher')
  METHOD = 2;
elseif strcmp(P.sig_calculation_method, 'LLRT')
  METHOD = 3;
elseif strcmpi(P.sig_calculation_method, 'concatenation')
  if P.use_1D_concatenation
    METHOD = 4;
  else
    % concatenation with categories
    if P.use_sample_specific_mutation_rates
      METHOD = 5;
    else
      METHOD = 6;
      fprintf('Substituting "concatenation"->"newcat" in P.sig_calculation_method\n');
      METHOD = 7;
    end
  end
elseif strcmpi(P.sig_calculation_method, 'newcat')
  METHOD = 7;
else
  error('Unknown significance calculation method "%s"', P.sig_calculation_method);
end

switch METHOD

case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  PROJECTION METHOD
%
%     Analyzes mutational data taking into account the possibility that
%     multiple mutations in the same gene in the same sample may not have
%     been independent events.   -- added 2008-04-23
%
%  1. Projects mutational profile of each sample to a restricted subspace,
%     where at most one mutation is represented,
%     selected from the most unlikely category.
%  2. For each sample, calculates the probability of that sample projecting
%     to each point in the collapsed projection.
%  3. Calculates P value for the gene based on the probability of obtaining
%     a distribution of samples in the collapsed space at least as extreme
%     as the observed one.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('PROJECTION METHOD NEEDS TO BE DEBUGGED.  See cases of TP53 and AMELY in ovarian dataset');

fprintf('Using projection method.\n');

silence = false;

nps = length(P.patient_subset);

for g=P.first_gene_to_calculate:P.last_gene_to_calculate

  if ~mod(g,1000), fprintf('%d/%d ', g, M.ng); end

  gname = M.gene.name{g};
%  if ~strcmp(gname,'TP53') && ~strcmp(gname,'AMELY'), continue; end
%  keyboard;

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  mu = repmat(rate(g,:),nps,1);
  N = squeeze(N_work(g,:,:))';
  Pmut = 1-((1-mu).^N);
  [Pmut priority] = sort(Pmut,2,'descend');
  Pclear = 1-Pmut;

  % STEP 2
  % for each sample, compute probability that it would been of each degree.

  Pdeg = ones(nps,nmutcats+1);
  for d=0:nmutcats
    % has to be clear in all categories > d
    Pdeg(:,d+1) = prod(Pclear(:,d+1:end),2);
    % and nonclear in category d (if d>0)
    if d>0
      Pdeg(:,d+1) = Pdeg(:,d+1) .* Pmut(:,d);
    end
  end

  % STEP 3
  % determine actual degree of each sample

  degree = zeros(nps,1);
  for p = 1:nps
    for c=fliplr(priority(p,:))
      if n_work(g,c,p)>0
         degree(p) = c;
         break
      end
    end
  end

  fprintf('WARNING: n_used has been abolished--projection method no longer reports correct tallies\n');
%  for c=1:nmutcats
%    n_used(g,c) = sum(degree==c);
%  end

  if P.zero_top_genes
    % make exception for top two genes
    if ismember(gname, P.top_genes_to_zero)
      Prob(g) = 0;
      s_g(g) = Inf;
      continue;   % don't do iterations for these genes
    end
  end

  % STEP 4
  % compute P value for gene by semiexact convolution procedure

  Prob(g) = semiexactP(Pdeg,degree,P.semi_exact_procedure_bin_size);

  % STEP 4a
  % compute score s_g to report in table

  Pobs = 1;
  for i=1:nps
    Pobs = Pobs * Pdeg(i,degree(i)+1);
  end
  if Pobs>0
    s_g(g) = -log10(Pobs);
  else
    s_g(g) = Inf;
  end

  % DONE

  if ~silence
    fprintf('%d\t%s\t%d\n', g, gname, Prob(g));
    %keyboard
  end

end % next gene

%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of projection method

case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  FISHER method
%
%      Per-category Fisher's Exact Test followed by Fisher Method of combining p-values
%
%  added 2010-07-21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using Fisher method.\n');

silence = true;

N_work = sum(N_work,3);
n_work = sum(n_work,3);

for g=P.first_gene_to_calculate:P.last_gene_to_calculate

  if ~mod(g,1000), fprintf('%d/%d ', g, M.ng); end

  gname = M.gene.name{g};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % calculate Fisher's Exact Test for each category
  %
  %               non-mutated   mutated
  %   this gene       [a]         [b]
  %   all genes       [c]         [d]
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  n = n_work(g,1:nmutcats);
  N = N_work(g,1:nmutcats); 
  p = nan(nmutcats,1);

  for i=1:nmutcats
    a = N(i)-n(i);              b = n(i);
    c = rate_N(i)-rate_n(i);    d = rate_n(i);
    if (b/a < d/c)
       p(i) = 1;    % hypomutated: doesn't count toward significance
    else
       p(i) = fisher_exact_test(a,b,c,d);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2
  % combine p-values using Fisher Method
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  X2 = -2 * sum(log(p));
  Prob(g) = 1 - chi2cdf(X2,2*nmutcats);

  if ~silence
      fprintf('%d\t%s\t%d\n', g, gname, Prob(g));
  end

 end   % next gene

fprintf ('Done\n');      % all genes finished!
%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of Fisher method


case 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LLRT METHOD
%
%      Log likelihood ratio test
%
%  added 2008-05-13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using LLRT.\n');

chatter = 0;
silence = 1;

N_work = sum(N_work,3);
n_work = sum(n_work,3);

for g=P.first_gene_to_calculate:P.last_gene_to_calculate

  if ~mod(g,1000), fprintf('%d/%d ', g, M.ng); end

  gname = M.gene.name{g};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % calculate point LLRT of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  n = n_work(g,1:nmutcats);
  N = N_work(g,1:nmutcats);
  mu_best = n ./ N;
  mu_null = rate(g,:);
  p_ratio = ((mu_best.^n).*((1-mu_best).^(N-n))) ./ ...
            ((mu_null.^n).*((1-mu_null).^(N-n)));
  p_ratio(isinf(p_ratio)) = 1e140;    % (happens for TP53 in ovarian data)
  p_ratio(isnan(p_ratio)) = 1;        % (happens with N=0)
  LLRT_obs = sum(log(p_ratio));
  if LLRT_obs<0, error('%s: LLRT < 0!\n', gname); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2
  % calculate P value of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if P.zero_top_genes
    % make exception for top two genes
    if ismember(gname, P.top_genes_to_zero)
      Prob(g) = 0;
      continue;   % don't do iterations for these genes
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2a
  % prepare table of piecewise probabilities
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = P.convolution_maxmuts;
  P_piece = zeros(nmutcats,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
  LLRT_piece = zeros(nmutcats,MAXMUTS+1);
  nmax = zeros(nmutcats+1,1);
  abort_flag = false;
  for c=1:nmutcats
    if abort_flag, break; end
    N = N_work(g,c);
    mu_null = rate(g,c);
    for n=0:MAXMUTS+1
      if abort_flag, break; end
      mu_best = n/N;
      Nchoosen = exp(ln_nchoosek(N,n));
      p_best = Nchoosen*(mu_best^n)*((1-mu_best)^(N-n));
      p_null = Nchoosen*(mu_null^n)*((1-mu_null)^(N-n));
      p_ratio = p_best / p_null;
      LLRT = log(p_ratio);
      P_piece(c,n+1)=p_null;
      LLRT_piece(c,n+1) = LLRT;
      if LLRT >= LLRT_obs    % this individual score would alone "break the bank"
         n = n - 1;
         break
      elseif N==0    % if zero coverage, do not consider making any mutations
         break
      elseif n > MAXMUTS
         fprintf('Exceeded MAXMUTS=%d with gene %d/%d (%s): setting P=0\n',MAXMUTS,g,M.ng,gname);
         abort_flag = true;
         Prob(g) = 0;
         continue;
      end
    end
    nmax(c) = n;
  end

  if abort_flag, continue; end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2b (SEMI-EXACT PROCEDURE)
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set up histogram

  hist_min_score = 0;
  hist_max_score = LLRT_obs;
  hist_bin_size = P.semi_exact_procedure_bin_size;
  hist_num_bins = ceil(hist_max_score / hist_bin_size) + 1;   % bin1 = zero only
  %if isnan(hist_num_bins) || isinf(hist_num_bins), keyboard; end
  H = zeros(hist_num_bins,1);


  % initial condition: all probability is in first bin (LLRT=0)

  H(1) = 1;

  % sequential convolution

  for c=1:nmutcats
    oldH = H;
    H = zeros(hist_num_bins,1);
    for n=0:nmax(c)
       piecewise_prob = P_piece(c,n+1);
       if (piecewise_prob == 0), continue; end
       piecewise_LLRT = LLRT_piece(c,n+1);
       nonempty_bins = find(oldH);
       for bin=1:length(nonempty_bins)
           bin_no = nonempty_bins(bin);
           old_prob = oldH(bin_no);
           LLRT_so_far = hist_bin_size * (bin_no-1);      % bin1 = zero only
           new_LLRT = LLRT_so_far + piecewise_LLRT;
           if new_LLRT < 0, new_LLRT = 0; end
           if new_LLRT < LLRT_obs
              new_prob = old_prob * piecewise_prob;
              new_bin_no = ceil(new_LLRT / hist_bin_size) + 1;   % bin1 = zero only
              H(new_bin_no) = H(new_bin_no) + new_prob;
           end
       end
    end
  end

  Pbulk = sum(H);
  Prob(g) = 1 - Pbulk;

  if ~silence
      fprintf('%d\t%s\t%d\n', g, gname, Prob(g));
  end

 end   % next gene

fprintf ('Done\n');      % all genes finished!
%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of LLRT method



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CONCATENATION METHOD
%
%     Analyzes mutational data by summing together all mutations in a given
%     gene over all samples, as if the gene is one long stretch of DNA
%     concatenated from the different patients.
%
%     Does not consider whether the mutations may be all from one patient.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ONE-DIMENSIONAL PROCEDURE
%
%  does not distinguish between mutation classes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using 1D concatenation method.\n');

mu = rate_tot;
for g=P.first_gene_to_calculate:P.last_gene_to_calculate, if ~mod(g,1000), fprintf('%d/%d ',g,M.ng); end
    N = Ntot_cov(g,M.TOT);
    n = ntot_nonsilent(g,M.TOT);
    Prob(g) = 1-binocdf(n-1,N,mu);
    if Prob(g)>0
      s_g(g) = -log10(Prob(g));
    else
      s_g(g) = Inf;
    end
end, fprintf('\n');
%%%%%%%%%%%%%%%%%%%% end of 1D procedure

case 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  CONCATENATION-WITH-CATEGORIES
%  with SAMPLE-SPECIFIC MUTATION RATES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using concatenation-with-categories with sample-specific mutation rate.\n');

silence = 0;

fprintf('gene ');

for g=P.first_gene_to_calculate:P.last_gene_to_calculate

  if mod(g,1000)==0, fprintf('%d/%d ', g, M.ng); end

  gname = M.gene.name{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:nmutcats,:))==0
    Prob(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % calculate point probability of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Pobs = 1;
  for patno=1:M.np
    for c=1:nmutcats
      N  = N_work(g,c,patno);
      n  = n_work(g,c,patno);
      mu = rate(patno,c);
      Pi = exp(ln_nchoosek(N,n))*(mu^n)*((1-mu)^(N-n));
      if isinf(Pi)||isnan(Pi), Pi = binopdf(n,N,mu); end
      Pobs = Pobs * Pi;
    end
  end

  % impose minimum value
  cutoff = 1e-100;
  if Pobs<cutoff, Pobs = cutoff; end
  if Pobs>1, error('Pobs>1'); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2
  % calculate P value of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if P.zero_top_genes
    % make exception for top two genes
    if ismember(gname, P.top_genes_to_zero)
      Prob(g) = 0;
      continue;   % don't do iterations for these genes
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2a
  % prepare table of piecewise probabilities
  % by patient, mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = P.convolution_maxmuts;
  Ppiece = zeros(M.np,nmutcats,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
  nmax = zeros(M.np,nmutcats+1);
  abort_flag = false;
  for patno=1:M.np
   for c=1:nmutcats
    if abort_flag, break; end
    N = N_work(g,c,patno);
    mu = rate(patno,c);
    for n=0:MAXMUTS+1
      if abort_flag, break; end
%      p = exp(ln_nchoosek(N,n))*(mu^n)*((1-mu)^(N-n));
%      if isinf(p)||isinf(p), p = binopdf(n,N,mu); end
      p = poisspdf(n,N*mu);
      Ppiece(patno,c,n+1) = p;
      if p < Pobs    % this individual score would alone "break the bank"
         n = n - 1;
         break
      elseif N==0    % if zero coverage, do not consider making any mutations
         break
      elseif Pobs==0 % this can happen with mu=0 in a random-subset permutation
         break
      elseif n > MAXMUTS
         fprintf('Exceeded MAXMUTS=%d with gene %d/%d (%s) pat %d/%d (%s): setting P=0\n',...
             MAXMUTS,g,M.ng,gname,patno,M.np,M.patient.name{patno});
         Prob(g) = 0;
         abort_flag = true;
         continue;
      end
    end
    nmax(patno,c) = n;
   end
 end

 if abort_flag, continue; end

 if Pobs>0
   score_obs = -log10(Pobs);
   s_g(g) = score_obs;
 else
   %%% (should no longer happen, because we're using a cutoff of 1e-100
   %    keyboard
   fprintf('Pobs of zero in concatenation-with-categories for gene %s!', gname);
   s_g(g) = 0;
 end

 hist_bin_size = P.semi_exact_procedure_bin_size;
  try_smaller_bin_size = true;
 while(try_smaller_bin_size)
  try_smaller_bin_size = false;

  % set up histogram
  hist_min_score = 0;
  hist_max_score = score_obs;
%  hist_bin_size = P.semi_exact_procedure_bin_size;
  hist_num_bins = ceil(hist_max_score / hist_bin_size) + 1;   % bin1 = zero only

  if isnan(hist_num_bins) || isinf(hist_num_bins)
    fprintf('ERROR with hist_num_bins, while processing gene %s\n', gname);
    keyboard
  end

  H = zeros(hist_num_bins,1);

  % initial condition: all probability is in first bin (score=0, P=1)
   H(1) = 1;

  try
    % sequential convolution
    for patno=1:M.np
      for c=1:nmutcats
        oldH = H;
        H = zeros(hist_num_bins,1);
        for n=0:nmax(patno,c)
          piecewise_prob = Ppiece(patno,c,n+1);
          if (piecewise_prob == 0), continue; end
          piecewise_score = -log10(piecewise_prob);
          nonempty_bins = find(oldH);
          for bin=1:length(nonempty_bins)
            bin_no = nonempty_bins(bin);
            old_prob = oldH(bin_no);
            score_so_far = hist_bin_size * (bin_no-1);      % bin1 = zero only
            new_score = score_so_far +  piecewise_score;
            if new_score < score_obs
              new_prob = old_prob * piecewise_prob;
              new_bin_no = ceil(new_score / hist_bin_size) + 1;   % bin1 = zero only
              H(new_bin_no) = H(new_bin_no) + new_prob;
    end,end,end,end,end

    Pbulk = sum(H);
    Prob(g) = 1 - Pbulk;
    if ~silence
      fprintf('%d\t%s\t%d\n', g, gname, Prob(g));
    end
    
  catch me
    disp(me);
    disp(me.message);
    fprintf('Error in convolutions with %s\n',gname);
    Prob(g) = 1;
  end

  if     Prob(g)<0.000001
    if hist_bin_size > P.semi_exact_procedure_bin_size/128
      hist_bin_size =  P.semi_exact_procedure_bin_size/128;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  elseif Prob(g)<0.00001
    if hist_bin_size > P.semi_exact_procedure_bin_size/64
      hist_bin_size =  P.semi_exact_procedure_bin_size/64;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  elseif Prob(g)<0.0001
    if hist_bin_size > P.semi_exact_procedure_bin_size/32
      hist_bin_size =  P.semi_exact_procedure_bin_size/32;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  elseif Prob(g)<0.001
    if hist_bin_size > P.semi_exact_procedure_bin_size/16
      hist_bin_size =  P.semi_exact_procedure_bin_size/16;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  elseif Prob(g)<0.01
    if hist_bin_size > P.semi_exact_procedure_bin_size/8
      hist_bin_size =  P.semi_exact_procedure_bin_size/8;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  elseif Prob(g)<0.1
    if hist_bin_size > P.semi_exact_procedure_bin_size/4
      hist_bin_size =  P.semi_exact_procedure_bin_size/4;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  elseif Prob(g)<0.5
    if hist_bin_size > P.semi_exact_procedure_bin_size/2
      hist_bin_size =  P.semi_exact_procedure_bin_size/2;
      fprintf('Retrying %s with bin_size = %f\n',gname,hist_bin_size);
      try_smaller_bin_size = true;
    end
  end
 end  % while (try_smaller_bin_size)

end % next gene
%%%% end of sample-specific-BMR method

case 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  MULTI-DIMENSIONAL PROCEDURE
%
%  COMPUTE P VALUE FOR EACH GENE
%  distinguishing between mutation classes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using concatenation-with-categories method.\n');

chatter = 0;
silence = 1;

N_work = sum(N_work,3);
n_work = sum(n_work,3);

fprintf('gene ');

for g=P.first_gene_to_calculate:P.last_gene_to_calculate

  if mod(g,1000)==0, fprintf('%d/%d ', g, M.ng); end
 
  gname = M.gene.name{g};
%  if strcmp(gname,'TP53'), keyboard; end
%disp(gname)

  % save time if gene has no mutations
  if fullsum(n_work(g,1:nmutcats))==0
    Prob(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % calculate point probability of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Pobs = 1;
  for c=1:nmutcats
    N  = N_work(g,c);
    n  = n_work(g,c);
    mu = rate(g,c);
    Pi = exp(ln_nchoosek(N,n))*(mu^n)*((1-mu)^(N-n));
    if isinf(Pi)||isnan(Pi), Pi = binopdf(n,N,mu); end
    Pobs = Pobs * Pi;
  end

  % impose minimum value
  cutoff = 1e-100;
  if Pobs<cutoff, Pobs = cutoff; end
  if Pobs>1, error('Pobs>1'); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2
  % calculate P value of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if P.zero_top_genes
    % make exception for top two genes
    if ismember(gname, P.top_genes_to_zero)
      Prob(g) = 0;
      continue;   % don't do iterations for these genes
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 2a
  % prepare table of piecewise probabilities
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  MAXMUTS = P.convolution_maxmuts;
  Ppiece = zeros(nmutcats,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
  nmax = zeros(nmutcats+1,1);
  abort_flag = false;
  for c=1:nmutcats
    if abort_flag, break; end
    N = N_work(g,c);
    mu = rate(g,c);
    for n=0:MAXMUTS+1
      if abort_flag, break; end
      p = exp(ln_nchoosek(N,n))*(mu^n)*((1-mu)^(N-n));
      if isinf(p)||isnan(p), p = binopdf(n,N,mu); end
      Ppiece(c,n+1) = p;
      if p < Pobs    % this individual score would alone "break the bank"
         n = n - 1;
         break
      elseif N==0    % if zero coverage, do not consider making any mutations
         break
      elseif Pobs==0 % this can happen with mu=0 in a random-subset permutation
         break
      elseif n > MAXMUTS
         fprintf('Exceeded MAXMUTS=%d with gene %d/%d (%s): setting P=0\n',MAXMUTS,g,M.ng,gname);
         Prob(g) = 0;
         abort_flag = true;
         continue;
      end
    end
    nmax(c) = n;
  end

  if abort_flag, continue; end

  if P.use_semi_exact_procedure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3a (SEMI-EXACT PROCEDURE)
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if Pobs>0
    score_obs = -log10(Pobs);
    s_g(g) = score_obs;
  else
   %%% (should no longer happen, because we're using a cutoff of 1e-100 
%    keyboard
    fprintf('Pobs of zero in concatenation-with-categories for gene %s!', gname);
    s_g(g) = 0;
  end

  % set up histogram

  hist_min_score = 0;
  hist_max_score = score_obs;
  hist_bin_size = P.semi_exact_procedure_bin_size;
  hist_num_bins = ceil(hist_max_score / hist_bin_size) + 1;   % bin1 = zero only

  if isnan(hist_num_bins) || isinf(hist_num_bins)
    fprintf('ERROR with hist_num_bins, while processing gene %s\n', gname);
  end

  H = zeros(hist_num_bins,1);

  % initial condition: all probability is in first bin (score=0, P=1)

  H(1) = 1;

  % sequential convolution

  for c=1:nmutcats
    oldH = H;
    H = zeros(hist_num_bins,1);
    for n=0:nmax(c)
       piecewise_prob = Ppiece(c,n+1);
       if (piecewise_prob == 0), continue; end
       piecewise_score = -log10(piecewise_prob);
       nonempty_bins = find(oldH);
       for bin=1:length(nonempty_bins)
           bin_no = nonempty_bins(bin);
           old_prob = oldH(bin_no);
           score_so_far = hist_bin_size * (bin_no-1);      % bin1 = zero only
           new_score = score_so_far +  piecewise_score;
           if new_score < score_obs
              new_prob = old_prob * piecewise_prob;
              new_bin_no = ceil(new_score / hist_bin_size) + 1;   % bin1 = zero only
              H(new_bin_no) = H(new_bin_no) + new_prob;
  end,end,end,end

  Pbulk = sum(H);
  Prob(g) = 1 - Pbulk;

  if ~silence
      fprintf('%d\t%s\t%d\n', g, gname, Prob(g));
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%% end of semi-exact procedure

  else % use exact procedure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3b (EXACT PROCEDURE)
  % examine all possible combinations of mutation counts
  %   for each, calculate a point probability Pposs.
  %   if Pposs > Pobs, then add Pposs to Pbulk.
  %   at end, 1-Pbulk is the P value of the observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Pbulk = 0;

  n = zeros(nmutcats+1,1);       % n(mutcats+1) is used as "done" flag
  nmax(nmutcats+1)=1; 

  if ~silence
    fprintf('gene %s\n', gname);
    if chatter
      nmax'
      % keyboard
    end
  end

  % BEGIN ITERATIONS

  first_term_flag = 1;

  while(~n(nmutcats+1))

    %
    % calculate point probability Pposs of current possibility
    %

    Pposs = 1;
    LPposs = 0;
    for c=1:nmutcats
      Pposs = Pposs * Ppiece(c,n(c)+1);
      LPposs = LPposs + LPpiece(c,n(c)+1);
    end

    if ~silence && chatter
          [n';nmax']
    end

    if Pposs > Pobs

    %
    % if it's in the "mountain" of more likely events, add it to Pbulk
    %

        Pbulk = Pbulk + Pposs;

        if first_term_flag
           LP1 = LPposs;
           S = 1;
           first_term_flag = 0;
        else
           if LPposs > LP1
              fprintf('Li>L1!  Maximum is not at 0,0,0,0!\n');
           end
           S = S + exp(LPposs-LP1);
        end

        % and the next configuration will simply be the next in series

        i = 1;                    % target first category for increment

        if ~silence && chatter
           fprintf('YES:  Pobs = %d Pposs = %d Pbulk = %d\n', Pobs, Pposs, Pbulk);
           keyboard
        end

    else

    %
    % else if it's NOT a more likely event, then cease incrementing this category
    %

        if any(n)                 % find next mutation category to increment
           i = find(n,1);         %     usually first nonzero category
        else                      % except if all-zeros "breaks the bank",
           i = nmutcats;          %     in which case, trigger end of iterations.
        end

        n(i) = 0;                 % first nonzero category gets zeroed out
        i = i + 1;                % next category will be targeted for increment

        if ~silence && chatter
          fprintf('NO:  Pobs = %d Pposs = %d Pbulk = %d\n', Pobs, Pposs, Pbulk);
          % keyboard
        end
    
    end

    %
    % now make the increment, carrying over to next category when necessary
    %

    while n(i) >= nmax(i)
       n(i) = 0;
       i = i + 1;
    end
    n(i) = n(i) + 1;

    if ~silence && ~chatter      % (i.e. unless full-verbose mode is on)
      if i>2
         [n';nmax']
         fprintf('%d %s Pobs = %d Pposs = %d Pbulk = %d\n', g, gname, Pobs, Pposs, Pbulk);
      end
    end

  end   % next iteration

% ITERATIONS FINISHED

  Pbulk = min(Pbulk,1);

  Prob(g) = 1 - Pbulk;

  if first_term_flag               % no conditions qualified
     PbulkL = 0;
  else
     LPbulk = LP1 + log(S);
     PbulkL = exp(LPbulk);
     PbulkL = min(PbulkL,1);
  end
  PL(g) = 1 - PbulkL;

  if ~silence
      %  fprintf('%d %s Pobs = %d Pbulk = %d PbulkL = %d P = %d PL = %d\n', ...
      %           g, gname, Pobs, Pbulk, PbulkL, P(g), PL(g));
      fprintf('%d\t%s\t%d\t%d\n', g, gname, Prob(g), PL(g));
      if abs((Prob(g)-PL(g))/Prob(g)) > 0.0000001, keyboard; end
  end

 end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of exact procedure

 end   % next gene

fprintf ('Done\n');      % all genes finished!

case 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% newcat
%
% new concatenation-with-categories method 2010/12/19
%
% --> performs "landfilling" to avoid problems with hypomutated cases
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Using "newcat" (landfilling) method.\n');

N_work = sum(N_work,3);
n_work = sum(n_work,3);

fprintf('gene ');

for g=P.first_gene_to_calculate:P.last_gene_to_calculate
  if mod(g,1000)==0, fprintf('%d/%d ', g, M.ng); end
  gname = M.gene.name{g};

  % save time if gene has no mutations
  if fullsum(n_work(g,1:nmutcats))==0
    Prob(g) = 1;
    continue
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 1
  % prepare table of piecewise probabilities and scores
  % by mutation category and number of mutations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MAXMUTS = 200;
  if any(n_work(g,:)*2>MAXMUTS)
    old = MAXMUTS;
    MAXMUTS = max(n_work(g,:))*2;
    fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
  end

  assign_p_zero = false;
  while(true)   % (may need to re-try with increased MAXMUTS)
    if MAXMUTS > P.convolution_maxmuts_limit
      fprintf('P.convolution_maxmuts_limit exceeded: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    prob_dist = zeros(nmutcats,MAXMUTS+1);    % second index is n+1 (to accommodate n=0)
    score_dist = zeros(nmutcats,MAXMUTS+1);
 
    for c=1:nmutcats
      N = N_work(g,c);
      mu = rate(g,c);
%      prob_dist(c,:) = binopdf(0:MAXMUTS,N,mu);
      prob_dist(c,:) = poisspdf(0:MAXMUTS,N*mu);
      score_dist(c,:) = -log10(prob_dist(c,:));
      nexp = floor(N*mu);
      if nexp>MAXMUTS
        old = MAXMUTS;
        MAXMUTS = nexp * 2;
        fprintf('Increasing MAXMUTS from %d to %d to cover number of observed mutations in %s\n',old,MAXMUTS,gname);
        continue;
      end
      if nexp>=1   % enforce one-sided scores ("landfill")        
        score_dist(c,1:nexp) = score_dist(c,nexp+1);
      end
    end
%      max_score = 500;  % <----- may need to change this treatment: can cause problem with summed scores
%      score_dist(score_dist>max_score) = max_score;
  % NOW WE JUST LEAVE IT AS Inf, and everything works out OK

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2
    % calculate score of the observation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    score_obs = 0;
    for c=1:nmutcats
      n  = n_work(g,c);
      score_obs = score_obs + score_dist(c,n+1);
    end

    if isinf(score_obs)
      fprintf('infinite score: assigning p=0\n');
      assign_p_zero = true;
      break;
    end

    % compute nmax for each category:
    %  the number of mutations where score_dist meets/exceeds score_obs
    need_to_retry = false;
    nmax = nan(nmutcats,1);
    for c=1:nmutcats
      idx = find(score_dist(c,:)>=score_obs,1);
      if isempty(idx)
        old = MAXMUTS;
        MAXMUTS = MAXMUTS * 2;
        fprintf('Increasing MAXMUTS from %d to %d to cover score_obs of %s\n',old,MAXMUTS,gname);
        need_to_retry = true;
        break;
      end
      nmax(c) = idx-1;
    end
    
    if ~need_to_retry, break; end
  end

  if assign_p_zero
    Prob(g) = 0;
    continue; % next gene
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 3
  % sequential convolution of mutation categories
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % set up histogram
  hist_min_score = 0;
  hist_max_score = score_obs;
  hist_bin_size = P.semi_exact_procedure_bin_size;
  hist_num_bins = ceil(hist_max_score / hist_bin_size) + 1;   % bin1 = zero only
  if isnan(hist_num_bins) || isinf(hist_num_bins)
    fprintf('ERROR with hist_num_bins, while processing gene %s\n', gname);
    keyboard
  end

  % initial condition: all probability is in first bin (score=0, P=1)
  H = zeros(hist_num_bins,1);
  H(1) = 1;

  % sequential convolution
  for c=1:nmutcats
    oldH = H;
    H = zeros(hist_num_bins,1);
    for n=0:nmax(c)
       piecewise_prob = prob_dist(c,n+1);
       if (piecewise_prob == 0), continue; end
       piecewise_score = score_dist(c,n+1);
       nonempty_bins = find(oldH);
       for bin=1:length(nonempty_bins)
           bin_no = nonempty_bins(bin);
           old_prob = oldH(bin_no);
           score_so_far = hist_bin_size * (bin_no-1);      % bin1 = zero only
           new_score = score_so_far +  piecewise_score;
           if new_score < score_obs
              new_prob = old_prob * piecewise_prob;
              new_bin_no = ceil(new_score / hist_bin_size) + 1;   % bin1 = zero only
              H(new_bin_no) = H(new_bin_no) + new_prob;
  end,end,end,end

  Pbulk = sum(H);
  Prob(g) = 1 - Pbulk;

 end   % next gene

end % switch(METHOD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% write results into M.gene
M.gene.score_classic = s_g;
M.gene.pval_classic = Prob;

% fix negative and unreliable very-small values
M.gene.pval_classic_lessthan_flag = (M.gene.pval_classic <= P.pval_cutoff);
M.gene.pval_classic(M.gene.pval_classic_lessthan_flag) = P.pval_cutoff;

fprintf ('Done\n');
