function M = load_mutdata2(project,patient_set,gene_set)
%
% loads data for a project (either 'TCGA' or 'TSP') into a structure
% and returns the structure
%
% Mike Lawrence 2008-06-17
%

if ~exist('project','var')
  error('"project" must be set to TCGA or TSP.\n');
end

if strcmp(project, 'TCGA')
if isnumeric(gene_set) && gene_set==601    % phase 1 only ("frozen" for publication)
  mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080702/tcga_20080711.mut';
  cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080502/context_20080606.txt';
else   % phases 1 and 2
  mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080812/tcga_20080819b.mut';
  cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080812/TCGA_context_20080813.txt';
end
elseif strcmp(project, 'TSP')
  mut_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_20080731.mut';
  cov_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_context_20080602.txt';
else
  error('"project" must be set to TCGA or TSP.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  OTHER INPUT FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% other input files

if strcmp(project, 'TCGA')
  if ~exist('gene_set','var'), gene_set = 601; end
  if isnumeric(gene_set)
    gene_file = ['/xchip/tcga/gbm/analysis/lawrence/tcga/' num2str(gene_set) 'genes.txt'];
  else
    gene_file = ['/xchip/tcga/gbm/analysis/lawrence/tcga/' gene_set];
  end
  if ~exist('patient_set','var'), patient_set = 94; end
  if isnumeric(patient_set)
    patient_file = ['/xchip/tcga/gbm/analysis/lawrence/tcga/' num2str(patient_set) 'patients.txt'];
  else
    patient_file = ['/xchip/tcga/gbm/analysis/lawrence/tcga/' patient_set];
  end
  targ_file = '/xchip/tcga/gbm/analysis/lawrence/cov/all_targets.fixed_names.txt';
  alias_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt';
  genecomp_file = '/xchip/tcga/gbm/analysis/lawrence/cov/2bp_gene_composition.txt';
  breakdown_file = '/xchip/tcga/gbm/analysis/lawrence/cov/TCGA_breakdown.txt';
  treated_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/treated_samples.txt';
  hypermutated_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/hypermutated_samples.txt';
  secondary_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/secondary_samples.txt';
  subtypes_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/subtypes.txt';
  purity_ploidy_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/purity_ploidy.txt';

elseif strcmp(project, 'TSP')
  if ~exist('gene_set','var'), gene_set = 623; end
  if isnumeric(gene_set)
    gene_file = ['/xchip/tcga/gbm/analysis/lawrence/tsp/' num2str(gene_set) 'genes.txt'];
  else
    gene_file = ['/xchip/tcga/gbm/analysis/lawrence/tsp/' gene_set];
  end
  if ~exist('patient_set','var'), patient_set = 188; end
  if isnumeric(patient_set)
    patient_file = ['/xchip/tcga/gbm/analysis/lawrence/tsp/' num2str(patient_set) 'patients.txt'];
  else
    patient_file = ['/xchip/tcga/gbm/analysis/lawrence/tsp/' patient_set];
  end
  targ_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_targets_2bp.txt';
  alias_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/aliases.txt';
  genecomp_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_composition_2bp.txt';
  breakdown_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/tsp_3N_stats.txt';
  rollup_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/new_TSP_rollup_FIXED.txt';
  hypermutated_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/hypermutated.txt';
  tertile_file = '/xchip/tcga/gbm/analysis/lawrence/tsp/sample_tertiles.txt';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  LOAD FILES AND PROCESS PATIENT + GENE NAMES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Loading input files for project %s:\n%s\n%s\n%s\n%s\n',...
   project, mut_file, cov_file, patient_file, gene_file);

M.file.mut = mut_file;
M.file.cov = cov_file;
M.file.patient = patient_file;
M.file.gene = gene_file;
M.file.breakdown = breakdown_file;
M.file.targ = targ_file;
M.file.alias = alias_file;

M.patient = tab2struct(patient_file);
M.gene = tab2struct(gene_file);
M.mut = tab2struct(mut_file);

M.cov = read_table(cov_file,'%s%s%f%f%f%f%f%f%f%f',char(9),1,'whitespace','\b\r');
M.genecomp = read_table(genecomp_file,'%s%f%f%f%f%f%f%f%f',char(9),1,'whitespace','\b\r');
M.breakdown = read_table(breakdown_file, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', ...
   char(9),1,'whitespace','\b\r');

if strcmp(project, 'TCGA')
  M.targ = tab2struct(targ_file,'%s%s%s%f%f%s%f');
  M.targ = rename_field(M.targ, 'name', 'gene');
  M.targ = rename_field(M.targ, 'ref', 'chr');
  M.targ = rename_field(M.targ, 'stop', 'end');
  M.targ = rename_field(M.targ, 'type', 'center');
elseif strcmp(project, 'TSP')
  M.targ = tab2struct(targ_file,'%s%s%f%f%s');
  M.rollup = read_table(rollup_file,['%s' repmat('%f', 1, 188)],char(9),1,'whitespace','\b\r');
end

% remove spliceform identifiers (eg .2)

M.targ.gene = regexprep(M.targ.gene,'\..*','');
M.breakdown.dat{1} = regexprep(M.breakdown.dat{1},'\..*','');

% replace old gene names with new gene names

aliases = tab2struct(alias_file);
M.gene.name = apply_aliases(M.gene.name,aliases);
M.targ.gene = apply_aliases(M.targ.gene,aliases);
M.genecomp.dat{1} = apply_aliases(M.genecomp.dat{1},aliases);
M.breakdown.dat{1} = apply_aliases(M.breakdown.dat{1},aliases);
M.cov.dat{2} = apply_aliases(M.cov.dat{2},aliases);
M.mut.gene = apply_aliases(M.mut.gene,aliases);

if strcmp(project, 'TSP')
  M.rollup.gene = apply_aliases(M.rollup.dat{1},aliases);
  M.rollup.patient = M.rollup.headers{1}(2:end)';
  M.rollup.cov = cat(2,M.rollup.dat{2:end});
end

M.breakdown.gene = M.breakdown.dat{1};
M.breakdown.stats = cat(2,M.breakdown.dat{2:end});
M.breakdown = rmfield(M.breakdown, {'headers','dat','dlm'});

M.mut.chr = convert_chr(M.mut.chr);

% extract phase+center for each TCGA gene

if strcmp(project,'TCGA')
  tmp = M.targ;
  idx = listmap(M.gene.name,tmp.gene);
  idx(isnan(idx)) = length(tmp.center)+1;
  tmp.center = [tmp.center; {''}];
  tmp.phase = [tmp.phase; NaN];
  M.gene.center = tmp.center(idx);
  M.gene.phase = tmp.phase(idx);
end

% map patient and gene names to indices to M.patient and M.gene

M.mut.patient = listmap(M.mut.patient, M.patient.name);
M.cov.patient = listmap(M.cov.dat{1}, M.patient.name);
if strcmp(project,'TSP')
  M.rollup.patient = listmap(M.rollup.patient, M.patient.name);
end

M.mut.gene = listmap(M.mut.gene, M.gene.name);
M.cov.gene = listmap(M.cov.dat{2}, M.gene.name);
M.genecomp.gene = listmap(M.genecomp.dat{1}, M.gene.name);
M.targ.gene = listmap(M.targ.gene, M.gene.name);
M.breakdown.gene = listmap(M.breakdown.gene, M.gene.name);
if strcmp(project, 'TSP')
  M.rollup.gene = listmap(M.rollup.gene, M.gene.name);
end

% Make sure mut data is properly marked "Do_not_use" for patients not in the patient list file

tmp = find(isnan(M.mut.patient));
M.mut.filtered(tmp) = repmat({'Do_not_use'},length(tmp),1);

M.ng = length(M.gene.name);
M.np = length(M.patient.name);

% load additional information when available

if exist('treated_file','var')
  tmp = load_struct(treated_file);
  tr = listmap(tmp.name,M.patient.name);
  M.patient.treated = false(M.np,1);
  tr = tr(~isnan(tr));
  M.patient.treated(tr) = true;
end

if exist('hypermutated_file','var')
  tmp = load_struct(hypermutated_file);
  hm = listmap(tmp.name,M.patient.name);
  M.patient.hypermutated = false(M.np,1);
  hm = hm(~isnan(hm));
  M.patient.hypermutated(hm) = true;
end

if exist('secondary_file','var')
  tmp = load_struct(secondary_file);
  sec = listmap(tmp.name,M.patient.name);
  M.patient.secondary = false(M.np,1);
  sec = sec(~isnan(sec));
  M.patient.secondary(sec) = true;
end

if exist('tertile_file','var')
  tmp = load_struct(tertile_file,'%s%s%f');
  idx = listmap(M.patient.name,tmp.name);
  idx(isnan(idx)) = length(tmp.tertile)+1;
  tmp.tertile = [tmp.tertile; nan];
  M.patient.tertile = tmp.tertile(idx);
end

if exist('subtypes_file','var')
  tmp = load_struct(subtypes_file,'%s%f');
  tmp.p = regexprep(tmp.SampleID,'_..._..$','');
  tmp.p = regexprep(tmp.p,'_','-');
  idx = listmap(M.patient.name,tmp.p);
  idx(isnan(idx)) = length(tmp.Cluster)+1;
  tmp.Cluster = [tmp.Cluster; NaN];
  M.patient.cluster = tmp.Cluster(idx);
end

if exist('purity_ploidy_file','var')
  tmp = load_struct(purity_ploidy_file,'%s%f%f%f');
  tmp.p = regexprep(tmp.sample,'-...-...$','');
  idx = listmap(M.patient.name,tmp.p);
  idx(isnan(idx)) = slength(tmp)+1;
  tmp.path_est_purity = [tmp.path_est_purity; NaN];
  tmp.est_purity = [tmp.est_purity; NaN];
  tmp.est_ploidy = [tmp.est_ploidy; NaN];
  M.patient.path_est_purity = tmp.path_est_purity(idx);
  M.patient.est_purity = tmp.est_purity(idx);
  M.patient.est_ploidy = tmp.est_ploidy(idx);
end

