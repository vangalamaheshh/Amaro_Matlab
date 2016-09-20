function M = analyze_mutation_data(M,outstem,P)
% M = analyze_mutation_data(M,outstem,P)
%
% Mike Lawrence 2010

if ~ischar(outstem), error('usage: analyze_mutation_data(M,outstem,P)'); end
if nargin<2, error('Needs M and outstem'); end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'using_scatter_gather',false);

P = impose_default_value(P,'quick',true);
if P.quick
  fprintf('QUICK MODE engaged: omitting some output files and analyses\n');
  P = impose_default_value(P,'output_per_exon_coverage',false);
  P = impose_default_value(P,'output_per_exon_mutations',false);
  P = impose_default_value(P,'output_per_gene_coverage',false);
  P = impose_default_value(P,'output_per_gene_mutation_counts',false);
  P = impose_default_value(P,'output_sample_sig_gene_table',false);
  P = impose_default_value(P,'perform_mutsig2_analysis',false);
  P = impose_default_value(P,'analyze_power_per_gene',false);
  P = impose_default_value(P,'max_genes_to_consider',10); % (for gene-gene correlations)
  P = impose_default_value(P,'number_of_top_genes_to_generate_coverage_plots_for',0);
  P = impose_default_value(P,'analyze_freq_genesets',false);
  P = impose_default_value(P,'analyze_clustered_mutations',false);
end

P = impose_default_value(P,'noplots',false);
if P.noplots
  P = impose_default_value(P,'generate_coverage_plot_and_bargraphs',false);
end
if ~isfield(P,'build') && isfield(M,'build') && ~isempty(M.build)
  P.build = M.build;
end
P = impose_default_value(P,'build','*required*');   % (for COSMIC and clustering analyses)
P = impose_default_value(P,'impute_full_coverage',~isfield(M,'N_cov'));
P = impose_default_value(P,'output_final_mutation_list',true);
P = impose_default_value(P,'output_per_gene_coverage',true);
P = impose_default_value(P,'output_per_gene_mutation_counts',true);
P = impose_default_value(P,'output_sample_sig_gene_table',true);
P = impose_default_value(P,'output_per_exon_coverage',true);
P = impose_default_value(P,'output_per_exon_mutations',false);
P = impose_default_value(P,'generate_coverage_plot_and_bargraphs',true);
P = impose_default_value(P,'perform_mutsig2_analysis',true);
P = impose_default_value(P,'skip_directly_to_mutsig2_analysis',false);
P = impose_default_value(P,'analyze_gene_gene_correlations',true);
P = impose_default_value(P,'analyze_cosmic_mutations',true);
P = impose_default_value(P,'analyze_clustered_mutations',true);
P = impose_default_value(P,'analyze_genesets',true);
P = impose_default_value(P,'analyze_freq_genesets',false);
P = impose_default_value(P,'analyze_power_per_gene',true);
P = impose_default_value(P,'number_of_top_genes_to_generate_coverage_plots_for',50);
P = impose_default_value(P,'interactive',false);

if P.perform_mutsig2_analysis
  MUTSIG_VERSION = '2.0';
else
  if P.use_sample_specific_mutation_rates && P.analyze_silent_nonsilent_ratios
    MUTSIG_VERSION = '1.5';
  else
    MUTSIG_VERSION = '1.0';
  end
end

if isfield(P,'geneset_collection_file') && ~isempty(P.geneset_collection_file) && ...
  ~strcmp(P.geneset_collection_file,'none') && ~strcmp(P.geneset_collection_file,'null')
  demand_file(P.geneset_collection_file);
end

if P.number_of_top_genes_to_generate_coverage_plots_for>0 && ~P.analyze_power_per_gene
  fprintf('Gene coverage plots need info from per-gene power analysis: will perform per-gene power analysis.\n');
  P.analyze_power_per_gene = true;
end

% if multi-mode, generate collapsed version for use in coverage plot and certain other analyses
if iscell(M)
  check_multiM_agreement(M);
  tmp = M;
  fprintf('Multi-set mode: collapsing sets for display\n');
  M = collapse_multiM(M);
  M.multi = tmp; clear tmp;
  multimode = true;
else 
  multimode = false;
end

if P.interactive
  fprintf('[interactive] '); keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    SUMMARY ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (summary analyses do not yet take advantage of scatter-gather)
if (~P.using_scatter_gather || P.is_gather) && ~P.skip_directly_to_mutsig2_analysis

  % output version report textfile
  outname = [outstem '.MutSig_version.txt'];
  save_textfile(['MutSig v' MUTSIG_VERSION],outname);
  % output mutation filtering report
  if isfield(M,'report') && ischar(M.report)
     outname = [outstem '.mutation_filtering_report.txt'];
    save_textfile(M.report,outname);
    fprintf('Saved %s\n',outname);
  end

  % draw coverage plot and bargraphs
  if P.generate_coverage_plot_and_bargraphs
    generate_coverage_plot_and_bargraphs(M,outstem,P)
  else
    fprintf('Skipping generation of coverage plot and bargraphs\n');
  end

  % output number of patients
  save_textfile(sprintf('%d\n',M.np),[outstem '.num_patients.txt']);

  % output final MAF file (only the mutations used in the analysis)
  if P.output_final_mutation_list
    outname = [outstem '.final_analysis_set.maf'];
    if isfield(M,'use')
      save_struct(reorder_struct(M.mut,M.use),outname);
    else
      save_struct(M.mut,outname);
    end
  end

  % output textfile with categories (for NicoGram)
  outname = [outstem '.mutcategs.txt'];
  save_struct(M.categ,outname);

  % output per-gene coverage and mutation counts
  if P.output_per_gene_coverage && ~P.impute_full_coverage
    outname = [outstem '.per_gene.coverage.txt'];
    output_gene_sample_coverage_table(M,outname);
  end
  if P.output_per_gene_mutation_counts
    outname = [outstem '.per_gene.mutation_counts.txt'];
    output_gene_sample_mutation_counts_table(M,outname);
  end

  % calculate mutation breakdown
  P.mutation_breakdown_output = [outstem '.mutation_breakdown.txt'];
  display_mutation_breakdown(M,P);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    MAIN SIGNIFICANCE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.siggene_report_filename = [outstem '.sig_genes.txt'];
P.mutation_rates_report_filename =  [outstem '.mutation_rates.txt'];

if P.output_sample_sig_gene_table
  P.sample_siggene_table_filename = [outstem '.sample_sig_gene_table.txt'];
end

if P.perform_mutsig2_analysis
  P.mutsig2_imagedir = [outstem '.mutsig2.images'];
end


M = find_significant_genes(M,P);  % (will handle scatter-gather if being used)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ADDITIONAL ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if multimode
  fprintf('Multimode implementation ends here.\n');
  return
end

if P.skip_directly_to_mutsig2_analysis
  fprintf('P.skip_directly_to_mutsig2_analysis=true:  skipping rest of analysis\n');
  return
end

if (P.using_scatter_gather && P.is_scatter)
  fprintf('Remaining analyses will be handled during gather.\n');

else
  M2 = M;
  if isfield(M2,'exprcorr')
    fprintf('Omitting expression-based BMR correction for remaining analyses\n');
    M2 = rmfield(M2,'exprcorr');
  end

  if isempty(P.geneset_collection_file) || strcmp(P.geneset_collection_file,'null') ||...
        strcmp(P.geneset_collection_file,'none') || ~exist(P.geneset_collection_file,'file')

    % don't perform genesets analysis
  else

    % analyze significantly mutated genesets
    if P.analyze_genesets
      P.geneset_analysis_excludes_genes = [];
      P.siggenesets_report_filename = [outstem '.sig_genesets.txt'];
      find_significant_genesets(M2,P); % (knows to disable scatter-gather)
      
      P.geneset_analysis_excludes_genes = M2.gene.name(M2.gene.qval<=0.1);
      % dependency note: has to wait until sig_genes analysis has been gathered!
      P.siggenesets_report_filename = [outstem '.sig_genesets_2.txt'];
      find_significant_genesets(M2,P); % (knows to disable scatter-gather)
    end

    % analyze frequently mutated genesets
    if P.analyze_freq_genesets
      P.freqgenesets_report_filename = [outstem '.freq_genesets.txt'];
      find_frequently_mutated_genesets(M2,P)
    end

  end

  if isempty(P.cosmic_file) || strcmp(P.cosmic_file,'null') ||...
        strcmp(P.cosmic_file,'none') || ~exist(P.cosmic_file,'file')

    % don't perform COSMIC analysis
  else
    try
      % analyze COSMIC mutations
      if P.analyze_cosmic_mutations
        P.cosmic_report_filename = [outstem '.cosmic_sig_genes.txt'];
        P.cosmic_report2_filename = [outstem '.cosmic_sig_genes2.txt'];
        P.cosmic_mutations_outname = [outstem '.cosmic_mutations.txt'];
        perform_COSMIC_overlap_analysis(M2,P);
      end
    catch me
      fprintf('ERROR in perform_COSMIC_overlap_analysis:\n%s\n',me.message);
    end
  end

  % analyze clustered mutations
  if P.analyze_clustered_mutations
    P.clustered_mutations_output = [outstem '.clustered_muts.txt'];
    try
      find_clustered_mutations(M2,P);
    catch me
      fprintf('ERROR in find_clustered_mutations:\n%s\n',me.message);
    end
  end
  
  % analyze mutual exclusivity / correlation of significantly mutated genes
  if P.analyze_gene_gene_correlations
    P.gene_gene_correlations_report_filename = [outstem '.gene_gene_correlations.txt'];
    P.gene_gene_correlations_verbose_report = false;
    % dependency note: has to wait until sig_genes analysis has been gathered!
    analyze_gene_gene_correlations_for_sig_genes(M2,P);
  end
  
  % analyze power per gene
  if P.analyze_power_per_gene 
    P.power_per_gene_outstem = [outstem '.power_per_gene'];
    try
      calculate_power_per_gene(M2,P)
    catch me
      fprintf('ERROR in calculate_power_per_gene:\n%s\n',me.message);
    end
  end
  
  % create coverage plots
  if P.number_of_top_genes_to_generate_coverage_plots_for>0
    P.per_gene_coverage_plots_outdir = [outstem '.gene_coverage_plots'];
    % dependency note: reports sig_gene info
    % dependency note: reports power per gene
    try
      generate_gene_coverage_plots(M2,P);
    catch me
      fprintf('ERROR in generate_gene_coverage_plots:\n%s\n',me.message);
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('analyze_mutation_data finished.\n');




