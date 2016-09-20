function M = find_significant_genesets(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'using_scatter_gather',false);

P=impose_default_value(P,'geneset_collection_file',...
  '/xchip/cga1/annotation/db/genesets/gsea_canonical_pathway_plus_HMT67.txt');
P=impose_default_value(P,'geneset_exclude','CANCER|MELANOMA|LEUKEMIA|GLIOMA|SCLEROSIS|DISEASE|DIABETES|CARCINOMA');
P=impose_default_value(P,'geneset_analysis_excludes_genes',{});

P=impose_default_value(P,'include_geneset_muttally_in_table',true);
P=impose_default_value(P,'include_geneset_description_in_table',true);
P=impose_default_value(P,'siggenesets_report_filename',[]);
P=impose_default_value(P,'sample_siggenesets_table_filename',[]);
P=impose_default_value(P,'siggenesets_tempfile_stem','sig_genesets');
P.siggene_report_filename = P.siggenesets_report_filename;
P.sample_siggene_table_filename = P.sample_siggenesets_table_filename;
P.siggene_tempfile_stem = P.siggenesets_tempfile_stem;
P.mutation_rates_report_filename = [];
P.use_precalculated_mutrates = true;
P.perform_mutsig2_analysis = false;

if P.using_scatter_gather
  if ~P.is_gather
    fprintf('This should only be called by the gather job!  Returning.\n');
    return;
  else
    P.using_scatter_gather = false;
    % this hack is necessary because the "sig genesets omitting sig genes" analysis
    % can't be run until the initial "sig genes" list has been gathered and finalized.
  end
end

M = map_to_genesets(M,P);
M = find_significant_genes(M,P);

