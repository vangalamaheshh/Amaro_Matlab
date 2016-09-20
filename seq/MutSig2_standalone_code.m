function MutSig2_standalone_code(mutation_file,coverage_file,covariate_file,...
                                 categ_file,target_file,coverage_fwb,context_and_effect_fwb,context_and_effect_categs_file,conservation_fwb,output_file)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Data-loading code from MutSigCV_v1_4_01
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % first, ensure output_file directory exists
  [outpath] = fileparts(output_file);
  if ~isempty(outpath) && ~exist(outpath,'dir'), mkdir(outpath); end
  
  % load MUTATION FILE and make sure all required fields are present

  fprintf('Loading mutation_file...\n');
  M = load_struct(mutation_file);

  % GENE
  if isfield(M,'gene') && isfield(M,'Hugo_Symbol')
    fprintf('NOTE:  Both "gene" and "Hugo_Symbol" are present in mutation_file.  Using "gene".\n');
  elseif isfield(M,'gene')
    % OK
  elseif isfield(M,'Hugo_Symbol')
    M.gene = M.Hugo_Symbol;
  else
    error('mutation_file lacks "gene" or "Hugo_Symbol" column.');
  end
  
  % PATIENT
  if isfield(M,'patient') && isfield(M,'Tumor_Sample_Barcode')
    fprintf('NOTE:  Both "patient" and "Tumor_Sample_Barcode" are present in mutation_file.  Using "patient".\n');
  elseif isfield(M,'patient')
    % OK
  elseif isfield(M,'Tumor_Sample_Barcode')
    M.patient = M.Tumor_Sample_Barcode;
  else
    error('mutation_file lacks "patient" or "Tumor_Sample_Barcode" column.');
  end

  % additional required columns in MAF
  if isfield(M,'is_coding') || isfield(M,'is_silent')
    fprintf('NOTE:  This version now ignores "is_coding" and "is_silent".  Requires "effect" column.\n');
    M = rmfield_if_exist(M,{'is_coding','is_silent'});
  end
  demand_field(M,{'categ','effect','chr','start','newbase','ref_allele'});  % (need to add correspondences to Chromosome, Start_position, etc.)
  M.chr = convert_chr(M.chr); M = make_numeric(M,'start');
  M.is_indel = grepmi('indel',M.effect) | strcmp('-',M.ref_allele) | strcmp('-',M.newbase);
  M.is_noncoding = strcmp('noncoding',M.effect);
  if isfield(M,'type')
    M.is_indel(grepmi('ins|del|frame',M.type)) = true;
    M.is_noncoding(grepmi('intron|utr|igr|flank|noncod',M.type)) = true;
  end
  M.newbase_idx = listmap(regexprep(M.newbase,'^(.).*$','$1'), {'A', 'C', 'G', 'T'});
  M.newbase_idx(M.is_indel) = 5;  % indels get newbase_idx=5

  % EFFECT
  M.effect = regexprep(M.effect,'^flank.*','noncoding');
  midx = find(M.is_indel & ~M.is_noncoding); M.effect(midx) = repmat({'indel_coding'},length(midx),1);
  
  effect_names = {'noncoding','silent','nonsilent','null','indel_coding'};
  M.effect_idx = listmap(M.effect,effect_names);
  if any(isnan(M.effect_idx))
    error('in mutation_file, "effect" must be one of noncoding/silent/nonsilent/null/indel_coding');
  end

  % LOAD COVERAGE FILE and COVARIATE FILE
  fprintf('Loading coverage file...\n');
  C = load_struct_specify_string_cols(coverage_file,1:3); % gene effect categ are all strings
  G=[]; G.gene = unique(C.gene); ng=slength(G);
  fprintf('Loading covariate file...\n');
  V = load_struct_specify_string_cols(covariate_file,1);  % gene is string
  f = fieldnames(V); cvnames = f(2:end); nv = length(cvnames);
  gidx = listmap(G.gene,V.gene);
  for i=1:length(f)
    if strcmp(f{i},'gene'), continue; end
    G.(f{i}) = nansub(V.(f{i}),gidx);
  end
  
  % make sure coverage file has required fields
  if ~isfield(C,'gene'), error('no "gene" column in coverage_file'); end
  if ~isfield(C,'effect') && isfield(C,'zone'), C = rename_field(C,'zone','effect'); end
  if ~isfield(C,'effect'), error('no "effect" column in coverage_file'); end
  C.effect = regexprep(C.effect,'^flank.*','noncoding');
  if any(~ismember(unique(C.effect),{'noncoding','silent','nonsilent'}))
    error('in coverage_file, "effect" must be one of noncoding/silent/nonsilent');
  end
  if ~isfield(C,'categ'), error('no "categ" column in coverage_file'); end
  f = fieldnames(C); coverage_patient_names = f(4:end);
  
  % remove any genes that we don't have coverage for
  badgene = setdiff(M.gene,C.gene);
  if ~isempty(badgene)
    fprintf('NOTE:  %d/%d gene names could not be mapped to coverage information.  Excluding them.\n',length(badgene),length(unique(M.gene)));
    M = reorder_struct_exclude(M,ismember(M.gene,badgene));
  end

  % make sure categories in C are same as in M
  bad = find(~ismember(M.categ,C.categ));
  if ~isempty(bad)
    fprintf('NOTE:  %d/%d mutations were outside the category set.  Excluding them.\n',length(bad),slength(M));
    tmp = sum(M.is_indel(bad)&M.is_noncoding(bad));
    if tmp>0, fprintf('\t(%d of them are noncoding indels.)\n',tmp); end
    M = reorder_struct_exclude(M,bad);
  end
  if slength(M)==0, error('No mutations left!\n'); end

  % map categories
  [K.name tmp C.categ_idx] = unique_keepord(C.categ);
  M.categ_idx = listmap(M.categ,K.name);
  if ~isfield(M,'categ_ignoring_null_categ'), M.categ_ignoring_null_categ = M.categ; end
  M.categ_ignoring_null_categ_idx = listmap(M.categ_ignoring_null_categ,K.name);
  ncat = slength(K);

  % ALSO load into K the full category specification information
  tmp = load_struct(categ_file);
  K = mapinto(K,tmp,'name');
  demand_fields(K,{'left','right','from','change','type'});

  % make sure there is a null+indel category
  knum = str2double(K.name);
  if all(~isnan(knum))
    % all category names are numbers: assume null+indel is the last one
    % (this maintains compatibility with the LUSC test data files)
    [tmp null_categ] = max(knum);
  else
    null_categ = grepi('null|indel',K.name,1);
    if length(null_categ)==0, error('ERROR: no null/indel category.\n'); end
    if length(null_categ)>1, error('ERROR: multiple null/indel categories.\n'); end
    % make sure all coding indels are assigned to this category
    midx = find(M.effect_idx==5);
    M.categ(midx) = repmat(K.name(null_categ),length(midx),1);
    M.categ_ignoring_null_categ(midx) = repmat(K.name(null_categ),length(midx),1);
    M.categ_idx(midx) = null_categ;
    M.categ_ignoring_null_categ_idx(midx) = null_categ;
  end

  % make sure C is sorted by the same gene order as in G
  C.gene_idx = listmap(C.gene,G.gene);
  C = sort_struct(C,'gene_idx');
  
  % map genes
  M.gene_idx = listmap(M.gene,G.gene);
  
  % regularize the sample name in the mutation table
  namebefore = M.patient;
  M.patient = regexprep(M.patient,'-Tumor$','');
  if any(~strcmp(namebefore,M.patient)), fprintf('NOTE:  Trimming "-Tumor" from patient names.\n'); end
  namebefore = M.patient;
  M.patient = regexprep(M.patient,'-','_');
  if any(~strcmp(namebefore,M.patient)), fprintf('NOTE:  Converting "-" to "_" in patient names.\n'); end

  pat=[]; [pat.name tmp M.patient_idx] = unique(M.patient);
  pat.cov_idx = listmap(pat.name,coverage_patient_names);
  np = slength(pat);
  if np<2, error('MutSig is not applicable to single patients.\n'); end

  % is generic coverage data given?
  generic_column_name = 'coverage';
  if length(coverage_patient_names)>1 || (length(coverage_patient_names)==1 && ~strcmpi(coverage_patient_names{1},generic_column_name))
    if any(strcmp(coverage_patient_names,generic_column_name)), error('reserved name "%s" cannot appear in list of patient names',generic_column_name); end
    if length(coverage_patient_names)~=length(unique(coverage_patient_names)), error('patient names in coverage_file must be unique'); end
    % make sure all patients are accounted for in coverage file
    if any(isnan(pat.cov_idx)), error('some patients in mutation_file are not accounted for in coverage_file'); end
    generic_coverage_flag = false;
  else
    % we're using generic coverage
    pat.cov_idx(:) = 1;
    generic_coverage_flag = true;
  end
  
  % BUILD n and N tables
  fprintf('Building n and N tables...\n');
  
  midx = ismember(M.effect_idx,find(strcmp('silent',effect_names)));
  n_silent = hist3d_fast(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  midx = ismember(M.effect_idx,grepmi('nonsilent|null|indel_coding',effect_names));
  n_nonsilent = hist3d_fast(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  midx = ismember(M.effect_idx,find(strcmp('noncoding',effect_names)));
  n_noncoding = hist3d_fast(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  % (change to hist3d for release)

  N_silent = nan(ng,ncat,np);
  N_nonsilent = nan(ng,ncat,np);
  N_noncoding = nan(ng,ncat,np);
  
  for ci=1:ncat
    silent_idx = strcmpi(C.effect,'silent') & C.categ_idx==ci;
    nonsilent_idx = strcmpi(C.effect,'nonsilent') & C.categ_idx==ci;
    noncoding_idx = strcmpi(C.effect,'noncoding') & C.categ_idx==ci;
    for pi=1:np
      cpfld = coverage_patient_names{pat.cov_idx(pi)};
      N_silent(:,ci,pi) = C.(cpfld)(silent_idx);
      N_nonsilent(:,ci,pi) = C.(cpfld)(nonsilent_idx);
      N_noncoding(:,ci,pi) = C.(cpfld)(noncoding_idx);
    end
  end

  % MAKE SURE ALL NUMBERS ARE INTEGERS
  n_silent = round(n_silent);
  n_nonsilent = round(n_nonsilent);
  n_noncoding = round(n_noncoding);
  N_silent = round(N_silent);
  N_nonsilent = round(N_nonsilent);
  N_noncoding = round(N_noncoding);

  % REMOVE MUTATIONS IN BINS WITH EXTREMELY LOW COVERAGE
  n_silent(n_silent>N_silent) = 0;
  n_nonsilent(n_nonsilent>N_nonsilent) = 0;
  n_noncoding(n_noncoding>N_noncoding) = 0;
  
  % SANITY CHECKS ON TOTALS
  tot_n_nonsilent = fullsum(n_nonsilent);
  tot_N_nonsilent = fullsum(N_nonsilent);
  tot_n_silent = fullsum(n_silent);
  tot_N_silent = fullsum(N_silent);
  tot_n_noncoding = fullsum(n_noncoding);
  tot_N_noncoding = fullsum(N_noncoding);
  tot_rate_nonsilent = tot_n_nonsilent/tot_N_nonsilent;
  tot_rate_silent = tot_n_silent/tot_N_silent;
  tot_rate_noncoding = tot_n_noncoding/tot_N_noncoding;
  tot_rate_coding = (tot_n_nonsilent+tot_n_silent)/(tot_N_nonsilent+tot_N_silent);
  
  min_tot_n_nonsilent = 50;
  min_tot_n_silent = 50;
  min_tot_n_noncoding = 50;
  min_rate_nonsilent = 1e-9;
  max_rate_nonsilent = 1e-3;
  min_rate_silent = 1e-9;
  max_rate_silent = 1e-3;
  min_rate_noncoding = 1e-9;
  max_rate_noncoding = 1e-3;
  max_abs_log2_difference_nonsilent_silent = 1.0;
  max_abs_log2_difference_noncoding_coding = 1.0;
  
  % see if silent and nonsilent are OK: if not, give warning
  if tot_n_nonsilent<min_tot_n_nonsilent || tot_n_silent<min_tot_n_silent, error('not enough mutations to analyze'); end
  if tot_rate_nonsilent<min_rate_nonsilent || tot_rate_nonsilent>max_rate_nonsilent, error('nonsilent mutation rate out of range'); end
  if tot_rate_silent<min_rate_silent || tot_rate_silent>max_rate_silent, error('silent mutation rate out of range'); end
  abs_log2_difference_nonsilent_silent = abs(log2(tot_rate_nonsilent/tot_rate_silent));
  if abs_log2_difference_nonsilent_silent>max_abs_log2_difference_nonsilent_silent, fprintf('Warning: silent and nonsilent rates are too different.\n'); end
  
  % see if noncoding is OK: if not, give warning and zero it all out
  ok = false;
  if tot_n_noncoding==0
    fprintf('NOTE:  no noncoding mutations.\n');
  else
    if tot_n_noncoding<min_tot_n_noncoding
      fprintf('WARNING:  not enough noncoding mutations to analyze\n');
    else
      if tot_rate_noncoding<min_rate_noncoding || tot_rate_noncoding>max_rate_noncoding
        fprintf('WARNING:  noncoding mutation rate out of range\n');
      else
        abs_log2_difference_noncoding_coding = abs(log2(tot_rate_noncoding/tot_rate_coding));
        if abs_log2_difference_noncoding_coding>max_abs_log2_difference_noncoding_coding
          fprintf('WARNING:  coding and noncoding rates are too different\n');
        else
          ok = true;
  end,end,end,end
  if ~ok
    fprintf('Zeroing out all noncoding mutations and coverage for the rest of the calculation.\n');
    n_noncoding(:) = 0;
    N_noncoding(:) = 0;
  end
  
  % add total columns
  n_silent(:,end+1,:) = sum(n_silent,2);
  n_nonsilent(:,end+1,:) = sum(n_nonsilent,2);
  n_noncoding(:,end+1,:) = sum(n_noncoding,2);
  N_silent(:,end+1,:) = N_silent(:,null_categ,:);          % copy total coverage from null+indel coverage
  N_nonsilent(:,end+1,:) = N_nonsilent(:,null_categ,:);
  N_noncoding(:,end+1,:) = N_noncoding(:,null_categ,:);
  
  % total across patients, save in G
  G.N_nonsilent = sum(N_nonsilent(:,end,:),3);
  G.N_silent = sum(N_silent(:,end,:),3);
  G.N_noncoding = sum(N_noncoding(:,end,:),3);
  G.n_nonsilent = sum(n_nonsilent(:,end,:),3);
  G.n_silent = sum(n_silent(:,end,:),3);
  G.n_noncoding = sum(n_noncoding(:,end,:),3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  MutSig2-specific code
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % MutSig2 parameters
  mutsig2_radius_to_impute_coverage_around_mutations = 10;
  mutsig2_max_coverage_bins = 20;
  mutsig2_randseed = 6789;
  mutsig2_maxperm = 2e5;
  mutsig2_theta = 1;
  mutsig2_min_effect_size = 1.01;
  mutsig2_clustering_metric = 204;
  mutsig2_restrict_to_one_mutation_per_patient = false;

  % initialize random number generator
  fprintf('Permutations: initializing with randseed=%d\n',mutsig2_randseed);
  rand('twister',mutsig2_randseed);

  % load list of target regions
  T = load_struct(target_file);
  demand_fields(T,{'gene','chr','start','end'});
  T.chr = convert_chr(T.chr);
  T = make_numeric(T,{'start','end'});
  T = reorder_struct_exclude(T,isnan(T.chr)|isnan(T.start)|isnan(T.end));
  if any(T.start>T.end), error('target list has start>end'); end
  T = sort_struct(T,{'gene','chr','start'});
  T.len = T.end-T.start+1;
  if mean(T.len)>5000, error('Looks like target list is whole transcripts.  Need individual exons!'); end
  T.gene_idx = listmap(T.gene,G.gene);

  % open FWB tracks
  FWB_conservation = org.broadinstitute.cga.tools.seq.FixedWidthBinary(conservation_fwb);
  FWB_conservation.setNullVal(200);
  FWB_context_and_effect = org.broadinstitute.cga.tools.seq.FixedWidthBinary(context_and_effect_fwb);
  context_and_effect_categs = load_struct(context_and_effect_categs_file);
  if ~strcmpi(coverage_fwb,'IMPUTE_FULL_COVERAGE')
    FWB_coverage = org.broadinstitute.cga.tools.seq.FixedWidthBinary(coverage_fwb);
    FWB_coverage.setNullVal(0);
  end

  % analyze collapsed categories vs. full-spectrum context_and_effect categories
  Q = collapse_context_and_effect_categories(context_and_effect_categs,K);

  % choose which mutations will be included in the permutations
  M.include_in_permutations = (M.effect_idx>=3 & M.effect_idx<=5 & M.categ_ignoring_null_categ_idx>=1 & M.categ_ignoring_null_categ_idx<=slength(K));

  % field for "position along track" that will be filled out during procedure
  M.trackpos = nan(slength(M),1);

  %%%%%%%%%%%%%%%%%%%%%%%
  % PROCESS EACH GENE
  %%%%%%%%%%%%%%%%%%%%%%%

  G.nperm = nan(ng,1);
  G.pCL = nan(ng,1);
  G.pFN = nan(ng,1);
  G.pCLFN = nan(ng,1);

keyboard

  for g=1:ng, if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
    try

    % PREPARATION FOR PERMUTATIONS

    % find the target regions for this gene
    tidx = find(T.gene_idx==g); genelength = sum(T.len(tidx));
    if genelength==0, continue; end        % no targets found

    % find the mutations for this gene and map to targets
    midx = find(M.gene_idx==g & M.include_in_permutations);
    exonstart=1;
    for ti=1:length(tidx),t=tidx(ti);
      z = find(M.chr(midx)==T.chr(t) & M.start(midx)>=T.start(t) & M.start(midx)<=T.end(t));
      M.trackpos(midx(z)) = exonstart + M.start(midx(z))-T.start(t);
      exonstart=exonstart+T.len(t);
    end
    midx(isnan(M.trackpos(midx)))=[];   % remove mutations that didn't map to a target

    if mutsig2_restrict_to_one_mutation_per_patient   % (chosen randomly)
      midx = midx(randperm(length(midx))); [u ui uj] = unique(M.patient_idx(midx)); midx = midx(ui);
    end    
    nm = length(midx);
    if nm<2, continue; end  % we only do permutations on genes with >=2 mutations
    
    % read conservation, coverage, and context_and_effect for these regions
    conservation_track = double(FWB_conservation.get(T.chr(tidx),T.start(tidx),T.end(tidx)));
    conservation_track(conservation_track==200) = NaN;  % missing data
    context_and_effect_track = double(FWB_context_and_effect.get(T.chr(tidx),T.start(tidx),T.end(tidx)));
    if ~strcmpi(coverage_fwb,'IMPUTE_FULL_COVERAGE')
      coverage_track = double(FWB_coverage.get(T.chr(tidx),T.start(tidx),T.end(tidx)));
    else
      coverage_track = ones(genelength,1);
    end

    % impute at least median coverage in radius around each mutation
    medcov = median(coverage_track(coverage_track>0));
    for i=1:nm,m=midx(i);
      i1 = max(1,M.trackpos(m)-mutsig2_radius_to_impute_coverage_around_mutations);
      i2 = min(genelength,M.trackpos(m)+mutsig2_radius_to_impute_coverage_around_mutations);
      coverage_track(i1:i2) = max(coverage_track(i1:i2),medcov);
    end

    % simplify coverage track (reduce number of discrete coverage levels)
    maxcov = max(coverage_track);
    if maxcov > mutsig2_max_coverage_bins
      coverage_track_factor = maxcov / mutsig2_max_coverage_bins;
      coverage_track = round(coverage_track / coverage_track_factor); % maybe should be ceil?
      maxcov = max(coverage_track);
    end
    if maxcov==0, continue; end   % this gene has no coverage

    % enumerate throwable positions for each mutation flavor
    [uce tmp throw_flavor] = unique([M.categ_ignoring_null_categ_idx(midx) M.effect_idx(midx)],'rows');
    nflavors = size(uce,1); flavor_counts = histc(throw_flavor,1:nflavors);
    throwable = cell(nflavors,1);
    for fi=1:nflavors
      c = uce(fi,1); e = uce(fi,2);
      if e<5          % e==3/nonsilent or e==4/null: need to match effect in Q
        [context_and_effect_idx newbase_idx] = find(Q(:,c,:)==e);
      else            % e==5/indel_coding: doesn't need to match syn/mis/nons, just has to be a "consistent" mutation, i.e. not A->A
        [context_and_effect_idx newbase_idx] = find(Q(:,c,:)>0);
        newbase_idx(:)=5; % ->indel
        % TO DO:  i think we can simplify this to just be the whole territory, weighted by the coverage track
      end
      cn = [context_and_effect_idx newbase_idx];
      % make sure the classes of the observed positions are included (in case track has disagreements)
      mm = midx(throw_flavor==fi);
      cn_obs = [context_and_effect_track(M.trackpos(mm)) M.newbase_idx(mm)];
      obs_only = setdiff(cn_obs,cn,'rows');
      if ~isempty(obs_only)
        mmm = nan(size(obs_only,1),1);
        for j=1:size(obs_only,1)
          mmm(j) = mm(find(context_and_effect_track(M.trackpos(mm))==obs_only(j,1) & M.newbase_idx(mm)==obs_only(j,2),1));
        end
        fprintf('___obs_only___     Gene %s  Flavor %d   Category %s   Effect %s\n',G.gene{g},fi,K.name{c},effect_names{e});
        pr(M.chr(mmm),M.start(mmm),M.type(mmm),M.splicedist(mmm),obs_only(:,1),context_and_effect_categs.name(obs_only(:,1)),obs_only(:,2));
      end
      cn = unique([cn;cn_obs],'rows');
      % make list of throwable for this flavor
      nx = size(cn,1); x = cell(nx,1);
      for i=1:nx
        idx = find(context_and_effect_track==cn(i,1));
        y = cell(length(idx),1);
        for j=1:length(idx), y{j} = repmat([idx(j) cn(i,2)],coverage_track(idx(j)),1); end
        x{i} = cat(1,y{:});
      end
      throwable{fi} = cat(1,x{:});
    end
    nthrowable = cellfun('length',throwable);
    if sum(nthrowable)==0, continue; end

    % range of metrics and bins for joint distrib
    min_clust = 0; max_clust = 1; min_cons = min(conservation_track); max_cons = max(conservation_track);
    nbins = 100; binsize_clust = (max_clust-min_clust)/(nbins-1); binsize_cons = (max_cons-min_cons)/(nbins-1);

    % calculate metrics for the observed mutations
    obs_cons_unreduced = nanmean(conservation_track(M.trackpos(midx)));
    obs_cons = obs_cons_unreduced / mutsig2_min_effect_size;
    obs_cons_bin = 1+floor((obs_cons - min_cons)/binsize_cons);
    obs_clust_unreduced = new_clustering_statistic(M.trackpos(midx),genelength,mutsig2_clustering_metric,true);
    obs_clust = obs_clust_unreduced / mutsig2_min_effect_size;
    obs_clust_bin = 1+floor((obs_clust - min_clust)/binsize_clust);

    % PERMUTATIONS

    k_clust = 0; k_cons = 0; joint_hist = zeros(nbins,nbins);
    nperm = 0; first_check = 1000; check_every = 10000; finished = false; tt1=tic;
    while(~finished), nperm = nperm + 1;
      % randomly throw mutations
      thrown_muts = cell(nflavors,1);
      for j=1:nflavors, if nthrowable(j)>0, thrown_muts{j} = throwable{j}(ceil(nthrowable(j)*rand(flavor_counts(j),1)),:);end,end
      perm_mutpos = cat(1,thrown_muts{:});

      % calculate metrics for this permutation and increment histograms
      perm_cons = nanmean(conservation_track(perm_mutpos(:,1)));
      perm_clust = new_clustering_statistic(perm_mutpos(:,1),genelength,mutsig2_clustering_metric);
      if perm_cons>=obs_cons, k_cons=k_cons+1; end
      if perm_clust>=obs_clust, k_clust=k_clust+1; end
      bin_cons = 1+floor((perm_cons - min_cons)/binsize_cons);
      bin_clust = 1+floor((perm_clust - min_clust)/binsize_clust);
      joint_hist(bin_clust,bin_cons) = joint_hist(bin_clust,bin_cons) + 1;
      
      if nperm~=first_check && mod(nperm,check_every)>0 && nperm<mutsig2_maxperm, continue; end  % don't need to check p-values every tieration

      % calculate marginal and joint p-values
      [p_clust ci_ratio_clust] = calc_pval_and_ci_ratio(k_clust,nperm);
      [p_cons ci_ratio_cons] = calc_pval_and_ci_ratio(k_cons,nperm);
      landfilled_hist = landfill(joint_hist/nperm);
      bin_score = -log10(landfilled_hist);
      obs_score = bin_score(obs_clust_bin,obs_cons_bin);
      k_joint = sum(joint_hist(bin_score>=obs_score));
      [p_joint ci_ratio_joint] = calc_pval_and_ci_ratio(k_joint,nperm);
      max_ci_ratio = max([ci_ratio_joint,ci_ratio_cons,ci_ratio_clust]);
      finished = (max_ci_ratio<=mutsig2_theta) | (nperm>=mutsig2_maxperm);

      fprintf('%s %5d/%5d %-11s nmuts %-4d len %-5d nperm %-7d (%0.f perm/sec) CONS k %-5d p %-0.6f  CLUST k %-5d p %-0.6f  JOINT k %-5d p %-0.6f\n',...
              datestr(now),g,ng,G.gene{g},nm,genelength,nperm,nperm/toc(tt1),k_cons,p_cons,k_clust,p_clust,k_joint,p_joint);
      keyboard
    end  % next permutation

    % record results of permutations
    G.nperm(g) = nperm;
    G.pCL(g) = p_clust;
    G.pFN(g) = p_cons;
    G.pCLFN(g) = p_joint;

    catch me
      fprintf('ERROR with gene %s\n',G.gene{g});
      disp(me);
      disp(me.message);
    end
  end   % next gene

  % close track files
  FWB_conservation.close();
  FWB_context_and_effect.close();
  if ~strcmpi(coverage_fwb,'IMPUTE_FULL_COVERAGE')
    FWB_coverage.close();
  end

  % for testing, just use the MutSig2 p-value
  G.p = G.pCLFN;

  % FDR
  G.q = calc_fdr_value(G.p);

  G = sort_struct(G,'p');
  save_struct(G,output_file);

  fprintf('Done.  Wrote results to %s\n',output_file);

end




  








