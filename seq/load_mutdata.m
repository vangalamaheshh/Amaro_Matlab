% load_mutdata.m

if ~exist('project','var')
  error('"project" must be set to TCGA or TSP.\n');
end

if strcmp(project, 'TCGA')
  mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080616/tcga_20080624.mut';
  cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080502/context_20080606.txt';
  sig_out_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080616/tcgasig.txt';
elseif strcmp(project, 'TSP')
  mut_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_20080602.mut';
  cov_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_context_20080602.txt';
  sig_out_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tspsig1.txt';
else
  error('"project" must be set to TCGA or TSP.\n');
end

% temporarily commented out for compiling purposes
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp     % for read_sample_info_file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  OTHER INPUT FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% other input files

if strcmp(project, 'TCGA')
  msi_file = '/xchip/tcga/gbm/results/sample_info/Master_sample_info_20080516.txt';
  targ_file = '/xchip/tcga/gbm/analysis/lawrence/cov/all_targets.fixed_names.txt';
  alias_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt';
  genecomp_file = '/xchip/tcga/gbm/analysis/lawrence/cov/2bp_gene_composition.txt';
  breakdown_file = '/xchip/tcga/gbm/analysis/lawrence/cov/TCGA_breakdown.txt';
elseif strcmp(project, 'TSP')
  targ_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_targets_2bp.txt';
  genecomp_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_composition_2bp.txt';
  breakdown_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_3N_stats.txt';
  gene_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tspgenes.txt';
  alias_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/aliases.txt';
  patient_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/patient_list.txt';
  rollup_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/new_TSP_rollup_FIXED.txt';
  was250_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/was250.txt';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LOAD FILES AND PROCESS PATIENT + GENE NAMES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Processing input files for project %s\n', project);

if strcmp(project,'TCGA')
  mut = tab2struct(mut_file,'%s%f%f%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
elseif strcmp(project,'TSP')
  mut = tab2struct(mut_file,'%s%f%f%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s');
end

cov = read_table(cov_file,'%s%s%f%f%f%f%f%f%f%f',char(9),1,'whitespace','\b\r');
genecomp = read_table(genecomp_file,'%s%f%f%f%f%f%f%f%f',char(9),1,'whitespace','\b\r');
breakdown = read_table(breakdown_file, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', ...
   char(9),1,'whitespace','\b\r');

if strcmp(project, 'TCGA')
  targ = read_table(targ_file, '%s%s%s%d%d%s%d',char(9),1,'whitespace','\b\r');
elseif strcmp(project, 'TSP')
  targ = read_table(targ_file, '%s%s%f%f%s',char(9),1,'whitespace','\b\r');
  gene = tab2struct(gene_file);
  rollup = read_table(rollup_file,['%s' repmat('%f', 1, 188)],char(9),1,'whitespace','\b\r');
  was250 = tab2struct(was250_file);
end

% remove spliceform identifiers (eg .2)

targ.dat{1} = regexprep(targ.dat{1},'\..*','');
breakdown.dat{1} = regexprep(breakdown.dat{1},'\..*','');
if strcmp(project, 'TSP')
  gene.gene = regexprep(gene.gene,'\..*','');
end

% replace old names with new names

aliases = tab2struct(alias_file);
targ.dat{1} = apply_aliases(targ.dat{1},aliases);
genecomp.dat{1} = apply_aliases(genecomp.dat{1},aliases);
breakdown.dat{1} = apply_aliases(breakdown.dat{1},aliases);
cov.dat{2} = apply_aliases(cov.dat{2},aliases);
mut.gene = apply_aliases(mut.gene,aliases);

if strcmp(project, 'TSP')
  rollup.gene = apply_aliases(rollup.dat{1},aliases);
  gene.gene = apply_aliases(gene.gene,aliases);
  was250.name = apply_aliases(was250.name,aliases);
  rollup.patient = rollup.headers{1}(2:end)';
  rollup.cov = cat(2,rollup.dat{2:end});
end

breakdown.gene = breakdown.dat{1};
breakdown.stats = cat(2,breakdown.dat{2:end});
breakdown = rmfield(breakdown, {'headers','dat','dlm'});

mut.chr = convert_chr(mut.chr);

if strcmp(project, 'TCGA')   % load master sample info file, find which patients are to be used
  SI=read_sample_info_file(msi_file);
  idx=grep('ABI',{SI.data_type},1);
  useme_x = find(strcmp({SI(idx).useme}', 'x'));
  tmp = {SI(idx).participant_id}';
  useme_names = unique(tmp(useme_x));
end

if strcmp(project, 'TSP')    % gene/center list
  [gene.name gi gj] = unique(gene.gene);
  gene.center = gene.center(gi);
  gene = rmfield(gene, 'gene');
end

% make comprehensive lists of patients and genes

patient_name = unique([mut.patient; cov.dat{1}]);
if strcmp(project, 'TSP')
  tmp = tab2struct(patient_file);
  patient_name = unique([patient_name; tmp.patient_id; rollup.patient]);
elseif strcmp(project, 'TCGA')
  patient_name = unique([patient_name; useme_names]);
end
mut.patient = listmap(mut.patient, patient_name);
cov.patient = listmap(cov.dat{1}, patient_name);
if strcmp(project, 'TCGA')
  useme = listmap(useme_names, patient_name);
elseif strcmp(project,'TSP')
  rollup.patient = listmap(rollup.patient, patient_name);
end

gene_name = unique([mut.gene; genecomp.dat{1}; targ.dat{1}; breakdown.gene; cov.dat{2}]);
if strcmp(project, 'TSP')
  gene_name = unique([gene_name; rollup.gene]);
end
mut.gene = listmap(mut.gene, gene_name);
cov.gene = listmap(cov.dat{2}, gene_name);
genecomp.gene = listmap(genecomp.dat{1}, gene_name);
targ.gene = listmap(targ.dat{1}, gene_name);
breakdown.gene = listmap(breakdown.gene, gene_name);
if strcmp(project, 'TSP')
  rollup.gene = listmap(rollup.gene, gene_name);
end


