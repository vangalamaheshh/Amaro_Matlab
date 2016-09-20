function M = load_all_mutation_data(P)
% Mike Lawrence 2009-09-14

if ~exist('P','var'), P=[]; end
if isfield(P,'whole_exome_genelist_input_file')
  P.genelist_input_file=P.whole_exome_genelist_input_file;
end
P=impose_default_value(P,'genelist_input_file',[]);
P=impose_default_value(P,'combined_mutation_list_input_file',[]);
P=impose_default_value(P,'combined_patient_list_input_file',[]);
P=impose_default_value(P,'ABI_coverage_input_file',[]);
P=impose_default_value(P,'capture_coverage_input_file',[]);
P=impose_default_value(P,'wgs_coverage_input_file',[]);
P=impose_default_value(P,'replace_missing_coverage',false);
P=impose_default_value(P,'coverage_replace_cutoff',0.05);
P=impose_default_value(P,'impose_average_cutoff',0.8);
P=impose_default_value(P,'remove_impossible_mutations',false);

try

P.project='TCGA';
P.mut_file = P.combined_mutation_list_input_file;
P.cov_file = P.ABI_coverage_input_file;
P.gene_file = P.genelist_input_file;
P.patient_file = P.combined_patient_list_input_file;

M = load_mutdata3(P);

P = impose_default_value(P,'use_all_mutations', true);
M = choose_mutations(M,P);
M = build_n_and_N_tables(M);
M = process_3N_breakdown_stats(M);

% replace_missing_coverage

if P.replace_missing_coverage, M = replace_missing_coverage(M,P); end

% collapse categories
M = collapse_mutation_categories(M);

% add capture coverage

if ~isempty(P.capture_coverage_input_file)
  fprintf('Adding capture sequencing coverage\n');
  tmp = load(P.capture_coverage_input_file);
  numcategs = 3;   % 4 in file, but only 3 are used: last one ("N") is always zero

  pidx = listmap(tmp.COV.patients,M.patient.name);
  pp = find(~isnan(pidx));
  gidx = listmap(tmp.COV.genes,M.gene.name);
  gg = find(~isnan(gidx));

  for p=1:length(pp), for c=1:numcategs
    M.N_cov(gidx(gg),c,pidx(pp(p))) = max(M.N_cov(gidx(gg),c,pidx(pp(p))), tmp.COV.cov(gg,pp(p),c));
  end,end
end

% add WGS coverage

if ~isempty(P.wgs_coverage_input_file)
  fprintf('Adding WGS sequencing coverage\n');
  tmp = load(P.wgs_coverage_input_file);
  numcategs = 3;   % 4 in file, but only 3 are used: last one ("N") is always zero

  pidx = listmap(tmp.COV.patients,M.patient.name);
  pp = find(~isnan(pidx));
  gidx = listmap(tmp.COV.genes,M.gene.name);
  gg = find(~isnan(gidx));

  for p=1:length(pp), for c=1:numcategs
    M.N_cov(gidx(gg),c,pidx(pp(p))) = max(M.N_cov(gidx(gg),c,pidx(pp(p))), tmp.COV.cov(gg,pp(p),c));
  end,end
end

% re-total
M.N_cov(:,M.TOT-M.NUM_INDEL_CLASSES:M.TOT,:) = ...
    repmat(sum(M.N_cov(:,1:M.TOT-M.NUM_INDEL_CLASSES-1,:),2),...
    [1 M.NUM_INDEL_CLASSES+1 1]);

% remove impossible mutations
if P.remove_impossible_mutations
  fprintf('Removing impossible mutations\n');
  M = remove_impossible_mutations(M,2);
end

% make "start" and "end" numeric
M.mut = make_numeric(M.mut,{'start','end'});

% add long genenames
M.gene.longname = get_longnames(M.gene.name);


catch me; excuse(me); end
