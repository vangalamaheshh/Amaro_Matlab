function OV = load_all_OV_mutation_data(P)
% Mike Lawrence 2009-07-02

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'combined_mutation_list_input_file','*required*');
P=impose_default_value(P,'combined_patient_list_input_file','*required*');
P=impose_default_value(P,'ABI_coverage_input_file','*required*');
P=impose_default_value(P,'capture_coverage_input_file','*required*');
P=impose_default_value(P,'whole_exome_genelist_input_file','*required*');

try

P.project='TCGA';
P.mut_file = P.combined_mutation_list_input_file;
P.cov_file = P.ABI_coverage_input_file;
P.gene_file = P.whole_exome_genelist_input_file;
P.patient_file = P.combined_patient_list_input_file;

OV = load_mutdata3(P);

P = impose_default_value(P,'use_all_mutations', true);
OV = choose_mutations(OV,P);
OV = build_n_and_N_tables(OV);
OV = process_3N_breakdown_stats(OV);

% replace_missing_coverage

P = impose_default_value(P,'coverage_replace_cutoff',0.05);
P = impose_default_value(P,'impose_average_cutoff',0.8);
OV = replace_missing_coverage(OV,P);
OV = collapse_mutation_categories(OV);

fprintf('Pruning ABI coverage to phase1 genes in the 24-patient set\n');
p24 = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/ova/24patients.txt');
pidx = find(~ismember(OV.patient.name,p24.name));
OV.N_cov(:,:,pidx) = 0;     % patients outside the ABI24 have zero ABI coverage

g601 = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/601genes.txt');
gidx = find(~ismember(OV.gene.name,g601.name));
OV.N_cov(gidx,:,:) = 0;     % genes outside phase1 have zero ABI coverage

% add NGS coverage

fprintf('Adding Next-Gen sequencing coverage\n');
load(P.capture_coverage_input_file)    %  NG
numcategs = 3;

pidx = listmap(NG.patients,regexprep(OV.patient.name,'^TCGA-\d\d-',''));
pp = find(~isnan(pidx));
gidx = listmap(NG.genes,OV.gene.name);
gg = find(~isnan(gidx));

for p=1:length(pp), for c=1:numcategs
  OV.N_cov(gidx(gg),c,pidx(pp(p))) = max(OV.N_cov(gidx(gg),c,pidx(pp(p))), NG.cov(gg,pp(p),c));
end,end
% re-total
OV.N_cov(:,OV.TOT-OV.NUM_INDEL_CLASSES:OV.TOT,:) = ...
  repmat(sum(OV.N_cov(:,1:OV.TOT-OV.NUM_INDEL_CLASSES-1,:),2),...
  [1 OV.NUM_INDEL_CLASSES+1 1]);

% remove impossible mutations
OV = remove_impossible_mutations(OV,2);

% make "start" and "end" numeric
OV.mut = make_numeric(OV.mut,{'start','end'});

catch me; excuse(me); end
