function M = find_significant_genes(M,P)
% Mike Lawrence 2008-2012

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'sig_calculation_method','projection');
P=impose_default_value(P,'use_precalculated_mutrates',false);
P=impose_default_value(P,'use_sample_specific_mutation_rates',true);
P=impose_default_value(P,'make_sample_specific_mutation_rates_equal_to_global',false);
if strcmpi(P.sig_calculation_method,'projection')
  P=impose_default_value(P,'pval_cutoff',1e-15);
else
  P=impose_default_value(P,'pval_cutoff',1e-11);
end
P=impose_default_value(P,'apply_exprcorr',isfield(M,'exprcorr'));
P=impose_default_value(P,'using_scatter_gather',false);
P=impose_default_value(P,'impute_full_coverage',~isfield(M,'N_cov'));
P=impose_default_value(P,'analyze_silent_nonsilent_ratios',true);
P=impose_default_value(P,'perform_mutsig2_analysis',~isfield(M.gene,'genes'));
P=impose_default_value(P,'mutsig2_imagedir',[]);  % default: will not write images
P=impose_default_value(P,'mutsig2_randseed',1234);  % will initialize once, then do all genes in sequence
P=impose_default_value(P,'skip_directly_to_mutsig2_analysis',false);
P=impose_default_value(P,'silent_combine_method','estimate_nonsilent_passengers');
P=impose_default_value(P,'mutation_rates_report_filename',[]);
P=impose_default_value(P,'siggene_report_filename',[]);
P=impose_default_value(P,'sample_siggene_table_filename',[]);
P=impose_default_value(P,'siggene_tempfile_stem','sig_genes');

multimode = isfield(M,'multi');
if multimode
  nsets = length(M.multi);
  check_multiM_agreement(M.multi);
end

if P.perform_mutsig2_analysis && ~P.impute_full_coverage
  demand_file(M.file.summed_cov_track);
end

if isfield(P,'genes_to_analyze') && ~isempty(P.genes_to_analyze)
  gta = P.genes_to_analyze;
  if ischar(gta), gta={gta}; end
  if islogical(gta), gta=find(gta); end
  if iscellstr(gta)
      gta = listmap(gta,M.gene.name);
  elseif isnumeric(gta)
    % OK
  else
    error('unknown format for P.genes_to_analyze');
  end
else
  gta = 1:M.ng;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    get basic statistics for each gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% total N and n and nsil
if ~P.skip_directly_to_mutsig2_analysis
  if ~P.impute_full_coverage
    Ntot_cov = sum(M.N_cov,3);
  else
    Ntot_cov = M.N_terr .* M.np;  
  end
  M.gene.N = Ntot_cov(:,M.TOT);
  nnon_tot = sum(M.n_nonsilent(:,1:M.TOT,:),3);
  M.gene.n = nnon_tot(:,M.TOT);
  if isfield(M,'n_silent')
    nsil_tot = sum(M.n_silent(:,1:M.TOT,:),3);
    M.gene.nsil = nsil_tot(:,M.TOT);
  end

  % npat, nsite
  M = count_sites_and_patients(M,gta);

  % mutation breakdown by (tumor type) or (category)
  if multimode
    for si=1:nsets
      fieldname = ['n_' genfieldname(M.multi{si}.name)];
      fieldcontents = sum(M.multi{si}.n_nonsilent(:,M.multi{si}.TOT,:),3);
      M.gene = setfield(M.gene,fieldname,fieldcontents);
    end
  elseif strcmpi(P.sig_calculation_method,'concatenation') || strcmpi(P.sig_calculation_method,'projection')
    if strcmpi(P.sig_calculation_method,'projection')
      fprintf('Projection mode warning: n subtotals in output table still reflect total mutation counts, not # mutated patients\n');
    end
    for c=1:M.TOT-1
      fieldname = ['n_' num2str(c)];
      fieldcontents = nnon_tot(:,c);
      M.gene = setfield(M.gene,fieldname,fieldcontents);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    calculate mutation rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~P.use_precalculated_mutrates || ~isfield(M,'mutrate')) && ~P.skip_directly_to_mutsig2_analysis
  [M mutrates] = estimate_mutation_rates(M,P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    prepare input for "classic" p-value calculation (working N, n, mu tables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BUILD WORKING n and N TABLES

if ~P.skip_directly_to_mutsig2_analysis
  if ~multimode
    nmutcats = M.TOT-1;
    if ~P.impute_full_coverage
      N_work = M.N_cov(:,1:nmutcats,:);
    else
      N_work = repmat(M.N_terr(:,1:nmutcats),[1,1,M.np]);
    end
    n_work = M.n_nonsilent(:,1:nmutcats,:);
  else % multimode
    nmutcats = 0;
    for i=1:nsets
      nmutcats = nmutcats + (M.multi{i}.TOT-1);
    end
    N_work = zeros(M.ng,nmutcats,M.np);
    n_work = zeros(M.ng,nmutcats,M.np);
    coffset = 0;
    poffset = 0;
    for i=1:nsets
      this_np = M.multi{i}.np;
      this_ncat = M.multi{i}.TOT-1;
      if ~P.impute_full_coverage
        N_work(:,coffset+[1:this_ncat],poffset+[1:this_np]) = M.multi{i}.N_cov(:,1:this_ncat,:);
      else
        N_work(:,coffset+[1:this_ncat],poffset+[1:this_np]) = repmat(M.multi{i}.N_terr(:,1:this_ncat),[1,1,this_np]);
      end
      n_work(:,coffset+[1:this_ncat],poffset+[1:this_np]) = M.multi{i}.n_nonsilent(:,1:this_ncat,:);
      coffset = coffset + this_ncat;
      poffset = poffset + this_np;
    end
  end

  % BUILD WORKING MUTATION RATE TABLE (patient x categ)
  
  if P.use_sample_specific_mutation_rates
    if multimode
      fprintf('2011-08-06: BUILD WORKING MUTATION RATE TABLE has been changed extensively\n');
      fprintf('            sample-specific multimode still needs to be coded\n');
      keyboard
    else
      if P.make_sample_specific_mutation_rates_equal_to_global
        fprintf('WARNING: setting all sample-specific mutation rates equal to the global rate\n');
        mu_work = repmat(sum(M.mutrate.ss.n,1)./sum(M.mutrate.ss.N,1),M.np,1);
      else
        mu_work = M.mutrate.ss.rate;
      end
    end
  else % global BMR
    if multimode
      fprintf('2011-08-06: BUILD WORKING MUTATION RATE TABLE has been changed extensively\n');
      fprintf('            Please tread carefully with multimode and make sure it works properly.\n');
      keyboard
      mu_work = cell(nsets,1);
      for i=1:nsets, mu_work{i} = M.multi{i}.mutrate.rate; end
      mu_work = cat(2,mu_work{:});
    else
      demand_field(M.mutrate,{'tot','rel'});
      mu_work = M.mutrate.tot.hat * M.mutrate.rel(1:nmutcats);
      mu_work = repmat(mu_work,M.np,1);  % expand to patients
    end
  end
  
  % expand mu_work to genes
  mu_work = repmat(shiftdim(mu_work',-1),M.ng,1);         %%% CHANGED 3/9/12: make sure it works
  
  % EXPRESSION-BASED BMR CORRECTION
  if isfield(M,'exprcorr') && P.apply_exprcorr
    fprintf('Applying expression-based gene-specific BMR correction\n');
    if size(M.exprcorr,1)~=M.ng, error('M.exprcorr should have one row per gene'); end
    if size(M.exprcorr,2)~=nmutcats
      if all(all(bsxfun(@minus,M.exprcorr(:,2:end),M.exprcorr(:,1))==0))
        % columns are all the same
        fprintf('Adjusting number of columns in supplied M.exprcorr\n');
        M.exprcorr = repmat(M.exprcorr(:,1),1,nmutcats);
      else
        % columns are not all the same
        error('M.exprcorr should have one column per category');
      end
    end
    mu_work = bsxfun(@times,mu_work,M.exprcorr);
  end
  
  % OTHER GENE-SPECIFIC BMR CORRECTION (typically coming from NS/S analysis)
  if isfield(M.mutrate,'per_gene_BMR_correction')
    fprintf('Applying final per-gene BMR correction (typically from NS/S analysis)\n');
    if size(M.mutrate.per_gene_BMR_correction,1)~=M.ng
      error('M.mutrate.per_gene_BMR_correction should have one row per gene')
    end
    if size(M.mutrate.per_gene_BMR_correction,2)~=nmutcats
      error('M.mutrate.per_gene_BMR_correction should have one column per category');
    end
    mu_work = bsxfun(@times,mu_work,M.mutrate.per_gene_BMR_correction);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    analyze silent/nonsilent ratios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~P.skip_directly_to_mutsig2_analysis
  if P.analyze_silent_nonsilent_ratios
    if (~P.impute_full_coverage && ... 
        (~isfield(M,'N_sil_cov') || ~isfield(M,'N_non_cov') || ~isfield(M,'n_nonsilent_ignoring_null_categ'))) || ...
          (P.impute_full_coverage && ...
           (~isfield(M,'N_sil_terr') || ~isfield(M,'N_non_terr') || ~isfield(M,'n_nonsilent_ignoring_null_categ')))
      fprintf('Cannot analyze silent/nonsilent ratios: required fields missing!\n');
      P.analyze_silent_nonsilent_ratios = false;
    end
  else
    fprintf('Not analyzing silent/nonsilent ratios.\n');
  end
  
  if P.analyze_silent_nonsilent_ratios
    
    M = analyze_silent_nonsilent_ratios(M,P); % (only to show p-value in table)
    
    %%%% VARIOUS WAYS OF USING NS/S INFO
    if strcmp(P.silent_combine_method,'sigmoid')
      % do nothing (will be handled after "classic" calculation)
    elseif strcmpi(P.silent_combine_method,'increase_prior')
      % examine S/NS ratios to see if we need to increase the per-gene BMRs
      M = examine_silent_nonsilent_ratios_as_evidence_against_prior(M,P);
      M.gene.per_gene_BMR_correction = M.mutrate.per_gene_BMR_correction;
    elseif strcmpi(P.silent_combine_method,'bounded_mle')
      % use MLE(nonsilent) bounded by CI(silent)
      P = impose_default_value(P,'silent_ci_to_use',0.90);
      M = compute_nonsilent_mle_bounded_by_silent_ci(M,P);
      M.gene.per_gene_BMR_correction = M.mutrate.per_gene_BMR_correction;
    elseif strcmpi(P.silent_combine_method,'estimate_nonsilent_passengers')
      % arbitrarily nominate and neutralize nonsilent passengers
      n_work = estimate_nonsilent_passengers_from_silent(n_work,M,P);
    else
      error('unknown P.silent_combine_method');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     calculate "CLASSIC" MutSig p-value:
%                 probability that all nonsilent mutations are due to background model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~P.skip_directly_to_mutsig2_analysis

  fprintf('\nRunning "classic" MutSig:\n');

  P.gene_names = M.gene.name;
  P.patient_names = M.patient.name;
  
  if ~P.using_scatter_gather
    
    [M.gene.pval M.gene.score] = calculate_significance(N_work,n_work,mu_work,P);
    
  else   % scatter-gather
    
    jobstem = [P.siggene_tempfile_stem '.classic'];
    jobcount = P.scatter_gather_jobcount;
    whichjob = P.scatter_gather_whichjob;
    [first_gene last_gene] = calculate_work_splits(jobcount, M.ng);
    
    if P.is_scatter
      
      P.first_gene_to_calculate = first_gene(whichjob);
      P.last_gene_to_calculate = last_gene(whichjob);
      fprintf('Scatter job %d: processing genes %d - %d\n',whichjob,P.first_gene_to_calculate,P.last_gene_to_calculate);
      genelistname = sprintf('%s.part%012d.genelist.txt',jobstem,whichjob);
      save_lines(M.gene.name(P.first_gene_to_calculate:P.last_gene_to_calculate),genelistname);
      [M.gene.pval M.gene.score] = calculate_significance(N_work,n_work,mu_work,P);
      tmp = reorder_struct(M.gene,P.first_gene_to_calculate:P.last_gene_to_calculate);
      outname = sprintf('%s.part%012d.mat',jobstem,whichjob);
      save(outname,'tmp');
      
    elseif P.is_gather
      
      fprintf('Gather job: collecting partial results\n');
      system(sprintf('mv scatter*/%s.part*.mat .',jobstem));
      tmp2 = cell(jobcount,1);
      for i=1:jobcount
        inname = sprintf('%s.part%012d.mat',jobstem,i);
        load(inname,'tmp');
        tmp2{i} = tmp;
      end
      M.gene = concat_structs(tmp2);
      system(sprintf('rm %s.part*.mat',jobstem));
      system(sprintf('rm scatter*/%s.part*.genelist.txt',jobstem));
      
    else
      error('invalid scatter-gather parameters');
    end
    
  end
  
  if ~P.using_scatter_gather || P.is_gather
    
    % fix negative and unreliable very-small values
    M.gene.pval_lessthan_flag = (M.gene.pval <= P.pval_cutoff);
    M.gene.pval(M.gene.pval_lessthan_flag) = P.pval_cutoff;
    
    if P.analyze_silent_nonsilent_ratios
      if strcmp(P.silent_combine_method,'sigmoid')
        M.gene = rename_field(M.gene,'pval','pval_classic');
        M.gene = rename_field(M.gene,'pval_lessthan_flag','pval_classic_lessthan_flag');
        M = combine_classic_and_ns_s_p_values(M,P);
      else
        % do nothing
      end
    end
  end

else

  % (if skipped classic analysis)
  M.gene.pval_classic = nan(slength(M.gene),1);
  M.gene.pval_classic_lessthan_flag = false(slength(M.gene),1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    MutSig 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if P.perform_mutsig2_analysis

  fprintf('\nRunning MutSig 2:\n');

  if ~P.using_scatter_gather

    R = call_mutsig_cluster(M,P);

  else   % scatter-gather
    
    jobstem = [P.siggene_tempfile_stem '.mutsig2'];
    jobcount = P.scatter_gather_jobcount;
    whichjob = P.scatter_gather_whichjob;
    [first_gene last_gene] = calculate_work_splits(jobcount, M.ng);

    if P.is_scatter

      fprintf('Scatter job %d: processing genes %d - %d\n',whichjob,first_gene(whichjob),last_gene(whichjob));
      P.genes_to_process = M.gene.name(first_gene(whichjob):last_gene(whichjob));
      P.mutsig2_imagedir = sprintf('mutsig2_images.part%012d',whichjob);
      R = call_mutsig_cluster(M,P);
      tmp = reorder_struct(R,first_gene(whichjob):last_gene(whichjob));
      outname = sprintf('%s.part%012d.mat',jobstem,whichjob);
      save(outname,'tmp');
      
    elseif P.is_gather

      fprintf('Gather job: collecting partial results\n');
      ensure_dir_exists(P.mutsig2_imagedir);
      system(sprintf('mv scatter*/mutsig2_images.part*/* %s',P.mutsig2_imagedir));
      system(sprintf('rm -rf scatter*/mutsig2_images.part*'));
      system(sprintf('mv scatter*/%s.part*.mat .',jobstem));
      tmp2 = cell(jobcount,1);
      for i=1:jobcount
        inname = sprintf('%s.part%012d.mat',jobstem,i);
        load(inname,'tmp');
        tmp2{i} = tmp;
      end
      R = concat_structs(tmp2);
      system(sprintf('rm %s.part*.mat',jobstem));
      
    else
      error('invalid scatter-gather parameters');
    end

  end

  if ~P.using_scatter_gather || P.is_gather

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%   INTEGRATE MUTSIG2 SCORES INTO FINAL P-VALUE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M.gene = mapinto(M.gene,R,'name','gene_name');

    if ~isfield(M.gene,'pval_classic') && isfield(M.gene,'pval')
      M.gene = rename_field(M.gene,'pval','pval_classic');
    end
    if ~isfield(M.gene,'pval_classic_lessthan_flag') && isfield(M.gene,'pval_lessthan_flag')
      M.gene = rename_field(M.gene,'pval_lessthan_flag','pval_classic_lessthan_flag');
    end
    M.gene.pval = fisher_combined_p([M.gene.pval_classic M.gene.p_joint]);
    idx = find(isnan(M.gene.p_joint));
    M.gene.pval(idx) = M.gene.pval_classic(idx);
    M.gene.pval_lessthan_flag = M.gene.pval_classic_lessthan_flag;
    idx = find(M.gene.p_joint<P.pval_cutoff);
    M.gene.pval(idx) = P.pval_cutoff;
    M.gene.pval_lessthan_flag(idx) = true;

  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    multiple-hypothesis correction and ranking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~P.using_scatter_gather || P.is_gather

  % calculate Benjamini-Hochberg FDR (q) values
  M.gene.qval = calc_fdr_value(M.gene.pval);

  % rank genes
  if strcmpi(P.sig_calculation_method,'projection') && isfield(M.gene,'npat')
    [tmp a] = sort(M.gene.npat, 'descend');
  elseif isfield(M.gene,'n') && isfield(M.gene,'N')
    [tmp a] = sort(M.gene.n./M.gene.N, 'descend');
  else
    a = 1:slength(M.gene);
  end
  [tmp b] = sort(M.gene.pval(a));
  M.gene.rank(a(b),1) = 1:M.ng;
  M.gene = order_fields_first(M.gene,{'rank'});

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    write reports
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~P.using_scatter_gather || P.is_gather

    % M.gene
    if isfield(P,'outdir')
      fname = [P.outdir '/results.mat'];
      gene=M.gene;
      fprintf('Saved M.gene to %s\n',fname);
      save(fname,'gene');
    end

    % mutation rates report
    if ~isempty(P.mutation_rates_report_filename) && exist('mutrates','var') && isstruct(mutrates)
      save_struct(mutrates,P.mutation_rates_report_filename)
    end

    % sample x significant gene table
    if ~isempty(P.sample_siggene_table_filename)
      output_sample_sig_gene_table(M,P);
    end

    % significantly mutated genes report
    if ~isempty(P.siggene_report_filename)
      output_sig_table(M,P);
    end
    
end
