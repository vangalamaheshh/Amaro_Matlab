function OV = load_OV(P)
% load_OV(P)
%
% P can be a parameter structure
% or the name of a dataset:
%   jamboree[default]
%

%default_dataset = 'jamboree';
default_dataset = '20091117';

if ~exist('P','var'), P = default_dataset; end

if exist('P','var') && ischar(P)
  tag = P;
  P=[];
else
  tag = [];
end

if ~exist('P','var'), P=[]; end
P.project='TCGA';

if ~isempty(tag), switch(tag)
case 'jamboree'
  P.combined_mutation_list_input_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20090511/ovarian_combined_20090519.mut';
  P.ABI_coverage_input_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081126/cov/tcga_ova_20081202.cov';
  P.need_to_standardize_patient_names = false;
  load_method = 1;
case 'jamboree-fixed-DNPs'
  P.combined_mutation_list_input_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20090511/ovarian_combined_20090519_fixed_DNPs.mut';
  P.ABI_coverage_input_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081126/cov/tcga_ova_20081202.cov';
  P.need_to_standardize_patient_names = false;
  load_method = 1;
case '20091020'
  %%%%% these files were generated using the following script:
  %%%%% /xchip/cga1/lawrence/ov/analysis/20091020/run.m
  P.combined_mutation_list_input_file = '/xchip/cga1/lawrence/ov/analysis/20091020/ov_combined_20091020.mut';
  P.ABI_coverage_input_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081126/cov/tcga_ova_20081202.cov';
  P.combined_patient_list_input_file = '/xchip/cga1/lawrence/ov/analysis/20091020/ov_full_patient_list.txt';
  P.cluster_assignments_file = [];
  P.capture_coverage_input_file = '/xchip/cga1/lawrence/ov/analysis/20091020/ov_cov_20091020_cap.mat';
  P.wgs_coverage_input_file = [];
  P.whole_exome_genelist_input_file = '/xchip/tcga_scratch/lawrence/capture/we_genelist.txt';
  P.need_to_standardize_patient_names = true;
 load_method = 2;
case '20091117'
  %%%%% these files were generated using the following script:
  %%%%% /xchip/cga1/lawrence/ov/analysis/20091113_3centers/run.m
  P.combined_mutation_list_input_file = '/xchip/cga1/lawrence/ov/analysis/20091111_3centers/3centers_filtered_by_e7.mut';
  P.ABI_coverage_input_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081126/cov/tcga_ova_20081202.cov';
  P.combined_patient_list_input_file = '/xchip/cga1/lawrence/ov/analysis/20091111_3centers/3centers_patient_list.txt';
  P.cluster_assignments_file = '/xchip/cga1/lawrence/ov/genepages/20091116/Top2500_k4_clusters.txt';
  P.capture_coverage_input_file = '/xchip/cga1/lawrence/ov/analysis/20091111_3centers/3centers_cov.mat';
  P.wgs_coverage_input_file = [];
  P.whole_exome_genelist_input_file = '/xchip/tcga_scratch/lawrence/capture/we_genelist.txt';
  P.need_to_standardize_patient_names = true;
  load_method = 2;
otherwise
  error('Unknown dataset %s',tag);
end, end


% LOAD ALIASES

aliases = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');

% LOAD MUTATION DATA

mut = [];
mut.mut = load_struct(P.combined_mutation_list_input_file);
mut.cov = load_struct(P.ABI_coverage_input_file);
mut.mut.gene = apply_aliases(mut.mut.gene,aliases);
mut.cov.gene = apply_aliases(mut.cov.gene,aliases);

% LOAD COPY-NUMBER DATA

% original file:
%  /xchip/tcga_scratch2/gsaksena/gistic_runs_20090512/call_genes_090512/standard/
%         standard_summarized_genes_spreadsheet.080630.txt
%  [linked to as "/xchip/cga1/lawrence/ov/genepages/cn.txt"]
%
% preprocessing:
%
% > cat cn.txt | sed 's/\bInconsistent\b/Inf/g' > cn_fixed1.txt
%
% 2009-10-20
% new data: up to batch 17 for all centers
%
% /xchip/tcga_scratch2/gsaksena/gistic_runs_2009-07-31-a/call_genes/standard/summarized_genes_spreadsheet.080630.txt
% [linked and preprocessed as above]
%
% [old files renamed to *_20090512.txt]
%
% 2009-11-17
%
% original file:   
%  /xchip/tcga/ovarian/analysis/gaddy/091115_mpcbs_gistic/b_9_15/call_genes/standard/latest_standard_run.genes.txt
%  [preprocessed as above--but apparently no instances of Inconsistent/Inf]

fprintf('Loading copy number data\n');
temp_file.cn = ['/xchip/tcga/ovarian/analysis/gaddy/091115_mpcbs_gistic/b_9_15/call_genes/standard/' ...
                'latest_standard_run.genes.txt'];
actual_file = '/xchip/cga1/lawrence/ov/genepages/20091116/cn_fixed1.txt';
nc = get_num_cols(actual_file);
cn = read_table(actual_file,['%s%s' repmat('%f',1,nc-2)],char(9),1,'whitespace','\b\r');
cn.gene.name = apply_aliases(cn.dat{1},aliases);
cn.gene.desc = cn.dat{2};
cn.patient.name = cn.headers{1}(3:end)';
cn.data = cat(2,cn.dat{3:end});
% cn.data(isinf(cn.data))=NaN;   % Inf->Nan


% LOAD EXPRESSION DATA

% original files:
% /xchip/tcga_scratch2/ov_jamboree_2009/exp/
%    TCGA_Affy_U133A_batch9-11-12-13_level3.gct.txt
%    TCGA_Affy_HuEx_batch9-11-12-13_level3.gct.txt
%
% linked to as:
% /xchip/cga1/lawrence/ov/genepages/
%    exprU.txt
%    exprH.txt
%
% preproecssing:
%
% (1) removed first two lines, saved as expr?_fixed1.txt
%
% (2) > cat expr?_fixed1.txt | sed 's/\bNA\b/NaN/g' > expr?_fixed2.txt
%
% 2009-11-17
%
% original file (from DCC): TCGA_Batch9-15_median.txt
%
% preprocessing:
% (1) cat TCGA_Batch9-15_median.txt | sed 's/\bNA\b/NaN/g' > expr_fixed.txt 
% 

fprintf('Loading expression data\n');

if 0
temp_file.expr1 = '/xchip/tcga_scratch2/ov_jamboree_2009/exp/TCGA_Affy_U133A_batch9-11-12-13_level3.gct.txt';
temp_file.expr2 = '/xchip/tcga_scratch2/ov_jamboree_2009/exp/TCGA_Affy_HuEx_batch9-11-12-13_level3.gct.txt';
expr1file = '/xchip/cga1/lawrence/ov/genepages/20090521/exprU_fixed2.txt';
expr2file = '/xchip/cga1/lawrence/ov/genepages/20090521/exprH_fixed2.txt';

nc = get_num_cols(expr1file);
expr1 = read_table(expr1file,['%s' repmat('%f',1,nc-1)],char(9),1,'whitespace','\b\r');
expr1.gene.name = apply_aliases(expr1.dat{1},aliases);
expr1.patient.code = expr1.headers{1}(2:end)';
idx = grep('^TCGA-\d\d-\d\d\d\d',expr1.patient.code,1);
expr1.patient.code = expr1.patient.code(idx);
expr1.patient.name = regexprep(expr1.patient.code,'(TCGA-\d\d-\d\d\d\d)-.*$','$1');
expr1.data = cat(2,expr1.dat{idx+1});

nc = get_num_cols(expr2file);
expr2 = read_table(expr2file,['%s' repmat('%f',1,nc-1)],char(9),1,'whitespace','\b\r');
expr2.gene.name = apply_aliases(expr2.dat{1},aliases);
expr2.patient.code = expr2.headers{1}(2:end)';
idx = grep('^TCGA-\d\d-\d\d\d\d',expr2.patient.code,1);
expr2.patient.code = expr2.patient.code(idx);
expr2.patient.name = regexprep(expr2.patient.code,'(TCGA-\d\d-\d\d\d\d)-.*$','$1');
expr2.data = cat(2,expr2.dat{idx+1});
end

temp_file.expr = 'DCC: TCGA_Batch9-15_median.txt';
exprfile = '/xchip/cga1/lawrence/ov/genepages/20091116/expr_fixed.txt';

nc = get_num_cols(exprfile);
expr = read_table(exprfile,['%s' repmat('%f',1,nc-1)],char(9),1,'whitespace','\b\r');
expr.gene.name = apply_aliases(expr.dat{1},aliases);
expr.patient.code = expr.headers{1}(2:end)';
idx = grep('^TCGA-\d\d-\d\d\d\d',expr.patient.code,1);
expr.patient.code = expr.patient.code(idx);
expr.patient.name = regexprep(expr.patient.code,'(TCGA-\d\d-\d\d\d\d)-.*$','$1');
expr.data = cat(2,expr.dat{idx+1});



% LOAD METHYLATION DATA

% original files:
%   /xchip/tcga_scratch2/ov_jamboree_2009/meth/
%      DNA_Methylation_TCGA_Ovarian_179Samples.txt   (data)
%      jhu-usc.edu_OV.HumanMethylation27.1.adf.txt   (probe glossary)
%
% linked to as:
%   /xchip/cga1/lawrence/ov/genepages/
%      methyl_data.txt
%      methyl_probes.txt
%
% preprocessing:
% > cat methyl_data.txt | sed 's/\bNA\b/NaN/g' > methyl_data_fixed1.txt
%

temp_file.methyl_data = '/xchip/tcga_scratch2/ov_jamboree_2009/meth/DNA_Methylation_TCGA_Ovarian_179Samples.txt';
temp_file.methyl_probes = '/xchip/tcga_scratch2/ov_jamboree_2009/meth/jhu-usc.edu_OV.HumanMethylation27.1.adf.txt';

datafile = '/xchip/cga1/lawrence/ov/genepages/20090521/methyl_data_fixed1.txt';
probesfile = '/xchip/cga1/lawrence/ov/genepages/20090521/methyl_probes.txt';

probes = load_struct(probesfile);

nc = get_num_cols(datafile);
methyl = read_table(datafile,['%s' repmat('%f',1,nc-1)],char(9),1,'whitespace','\b\r');
idx = listmap(methyl.dat{1},probes.IlmnID);
if any(isnan(idx))
  fprintf('Error in processing methylation data:\n  some probes map to unknown gene\n');
  keyboard
end
methyl.gene.name = apply_aliases(probes.SYMBOL(idx),aliases);
methyl.patient.code = methyl.headers{1}(2:end)';
idx = grep('^TCGA-\d\d-\d\d\d\d',methyl.patient.code,1);
methyl.patient.code = methyl.patient.code(idx);
methyl.patient.name = regexprep(methyl.patient.code,'(TCGA-\d\d-\d\d\d\d)-.*$','$1');
methyl.data = cat(2,methyl.dat{idx+1});

% for each gene, take average methylation across probes for that gene
[g gi gj] = unique(methyl.gene.name);
keep = true(size(methyl.data,1),1);
for i=1:length(g)
  idx = find(gj==i);
  methyl.data(idx(1),:) = mean(methyl.data(idx,:));
  keep(idx(2:end)) = false;
end
methyl.data = methyl.data(keep,:);
methyl.gene = reorder_struct(methyl.gene,keep);

%%%%%
%%%%% MAP SEQUENCING DATA INTO MASTER MATRIX
%%%%%

% unified patient list

all_patients = [];
if P.need_to_standardize_patient_names
  mut.mut.patient = regexprep(mut.mut.patient,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
  mut.cov.patient = regexprep(mut.cov.patient,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
  methyl.patient.name = regexprep(methyl.patient.name,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
%  expr1.patient.name = regexprep(expr1.patient.name,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
%  expr2.patient.name = regexprep(expr2.patient.name,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
  expr.patient.name = regexprep(expr.patient.name,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
  cn.patient.name = regexprep(cn.patient.name,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
end
%all_patients.name = [mut.mut.patient; mut.cov.patient;...
%  methyl.patient.name;expr1.patient.name;expr2.patient.name;cn.patient.name];

all_patients.name = [mut.mut.patient; mut.cov.patient;...
  methyl.patient.name;expr.patient.name;cn.patient.name];

if isfield(P,'combined_patient_list_input_file')
  tmp = load_struct(P.combined_patient_list_input_file);
  tmp.name = regexprep(tmp.name,'TCGA-\d\d-(\d\d\d\d)','OV-$1');
  all_patients.name = [all_patients.name;tmp.name];
end
all_patients.name = unique(all_patients.name);

all_genes = [];
%all_genes.name = [mut.mut.gene; mut.cov.gene;...
%  methyl.gene.name;expr1.gene.name;expr2.gene.name;cn.gene.name];
all_genes.name = [mut.mut.gene; mut.cov.gene;...
  methyl.gene.name;expr.gene.name;cn.gene.name];

all_genes.name = unique(all_genes.name);

if load_method==1

  % reload mutation data using all_genes and all_patients  rr = num2str(round(rand*1e9));
  P.gene_file = ['/tmp/OV_genepages_all_genes_' rr '_.txt'];
  P.patient_file = ['/tmp/OV_genepages_all_patients_' rr '_.txt'];
  save_struct(all_patients,P.patient_file);
  save_struct(all_genes,P.gene_file);

  OV = load_mutdata3(P);
  OV.file = merge_structs({temp_file,OV.file});
  OV.mut.start = str2double(OV.mut.start);
  OV.mut.end = str2double(OV.mut.end);

  OV.mut.filtered = repmat({'OK'},slength(OV.mut),1);
  OV.mut.type = regexprep(OV.mut.type,'Read-through','Nonsense');   % (temporary fix)

  % mutation classes (categories)
  OV.mut.class = upper(OV.mut.ref_allele);
  %OV.mut.class = regexprep('^(C|G)$','$1 \(other\)');
  idx = union(grep('In|Del',OV.mut.type,1),find(OV.mut.end-OV.mut.start>0));
  OV.mut.class(idx) = repmat({'Indel'},length(idx),1);

 % splice-site indels --> frameshift     (temporary fix)
  idx = find(strcmp(OV.mut.type,'Splice_site') & strcmp(OV.mut.class,'Indel'));
  OV.mut.type(idx) = repmat({'Frameshift'},length(idx),1);

  for i=1:slength(OV.mut)
    if strcmp(OV.mut.class{i},'C')
      base = upper(genome_region(OV.mut.chr(i),OV.mut.start(i)+1));
      if strcmp(base,'G'), OV.mut.class{i} = 'C (CpG)';
      else OV.mut.class{i} = 'C (other)'; end
    elseif strcmp(OV.mut.class{i},'G')
      base = upper(genome_region(OV.mut.chr(i),OV.mut.start(i)-1));
      if strcmp(base,'C'), OV.mut.class{i} = 'G (CpG)';
      else OV.mut.class{i} = 'G (other)'; end
    end
  end

  OV.mut.type = regexprep(OV.mut.type,'Missense \(DNP\)','Missense');   % temporary fix

  P = impose_default_value(P,'use_all_mutations', true);
  OV = choose_mutations(OV,P);
  OV = build_n_and_N_tables(OV);
  OV = process_3N_breakdown_stats(OV);

  % replace_missing_coverage

  P = impose_default_value(P,'coverage_replace_cutoff',0.05);
  P = impose_default_value(P,'impose_average_cutoff',0.8);
  OV = replace_missing_coverage(OV,P);
  OV = collapse_mutation_categories(OV);

  p24 = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/ova/24patients.txt');
  pidx = find(~ismember(OV.patient.name,p24.name));
  OV.N_cov(:,:,pidx) = 0;     % patients outside the ABI24 have zero ABI coverage

  g601 = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/601genes.txt');
  gidx = find(~ismember(OV.gene.name,g601.name));
  OV.N_cov(gidx,:,:) = 0;     % genes outside phase1 have zero ABI coverage

  if 0
    gcov=sum(OV.N_cov(:,OV.TOT,:),3); sum(gcov>0)   % 591
    pcov=sum(OV.N_cov(:,OV.TOT,:),1); sum(pcov>0)   % 24
  end

  % add NGS coverage

  fprintf('Adding Next-Gen sequencing coverage\n');
  ngs_cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20090511/ngs_cov.mat';
  load(ngs_cov_file)    %  NG
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

  if 0
    gcov=sum(OV.N_cov(:,OV.TOT,:),3); sum(gcov>0)   % 17315
    pcov=sum(OV.N_cov(:,OV.TOT,:),1); sum(pcov>0)   % 53
  end

elseif load_method==2

  OV = load_all_mutation_data(P);
  OV.file = merge_structs({temp_file,OV.file});
else
  error('Unknown load_method');
end


%%%
%%%   MAP OTHER DATA MODALITIES INTO MASTER MATRIX
%%%

fprintf('Mapping all data to master matrix\n');

OV.gene.desc = repmat({''},OV.ng,1);
OV.patient.have_seq = squeeze(sum(sum(OV.N_cov(:,OV.TOT,:),1),2))>0;
OV.gene.have_seq = squeeze(sum(sum(OV.N_cov(:,OV.TOT,:),3),2))>0;

% map CN data

gidx = listmap(OV.gene.name,cn.gene.name);
pidx = listmap(OV.patient.name,cn.patient.name);
cn.data(:,end+1)=NaN;
cn.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(cn.data,1);
pidx(isnan(pidx))=size(cn.data,2);
OV.cn = cn.data(gidx,pidx);
OV.patient.have_cn = ~all(isnan(OV.cn),1)';
OV.gene.have_cn = ~all(isnan(OV.cn),2);

cn.gene.desc = [cn.gene.desc; {''}];
OV.gene.desc = cn.gene.desc(gidx);

% map EXPR data

if 0
gidx = listmap(OV.gene.name,expr1.gene.name);
pidx = listmap(OV.patient.name,expr1.patient.name);
expr1.data(:,end+1)=NaN;
expr1.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(expr1.data,1);
pidx(isnan(pidx))=size(expr1.data,2);
OV.expr1 = expr1.data(gidx,pidx);
OV.patient.have_expr1 = ~all(isnan(OV.expr1),1)';
OV.gene.have_expr1 = ~all(isnan(OV.expr1),2);

gidx = listmap(OV.gene.name,expr2.gene.name);
pidx = listmap(OV.patient.name,expr2.patient.name);
expr2.data(:,end+1)=NaN;
expr2.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(expr2.data,1);
pidx(isnan(pidx))=size(expr2.data,2);
OV.expr2 = expr2.data(gidx,pidx);
OV.patient.have_expr2 = ~all(isnan(OV.expr2),1)';
OV.gene.have_expr2 = ~all(isnan(OV.expr2),2);

% for now, only use expr2 (HuEx)
OV.expr = OV.expr2;
OV.patient.have_expr = OV.patient.have_expr2;
OV.gene.have_expr = OV.gene.have_expr2;
OV.file.expr = OV.file.expr2;
end

gidx = listmap(OV.gene.name,expr.gene.name);
pidx = listmap(OV.patient.name,expr.patient.name);
expr.data(:,end+1)=NaN;
expr.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(expr.data,1);
pidx(isnan(pidx))=size(expr.data,2);
OV.expr = expr.data(gidx,pidx);
OV.patient.have_expr = ~all(isnan(OV.expr),1)';
OV.gene.have_expr = ~all(isnan(OV.expr),2);


% map METHYL data

gidx = listmap(OV.gene.name,methyl.gene.name);
pidx = listmap(OV.patient.name,methyl.patient.name);
methyl.data(:,end+1)=NaN;
methyl.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(methyl.data,1);
pidx(isnan(pidx))=size(methyl.data,2);
OV.methyl = methyl.data(gidx,pidx);
OV.patient.have_methyl = ~all(isnan(OV.methyl),1)';
OV.gene.have_methyl = ~all(isnan(OV.methyl),2);


%%%
%%%  GENE METADATA
%%%


% PARSE GENE LOCATIONS

fprintf('Parsing gene locations\n');

Cyto = load_struct('/xchip/tcga/gbm/analysis/lawrence/genome/hg18/cytoBand.txt','%s%f%f%s%s');

OV.gene.chr = cell(OV.ng,1);
OV.gene.band = cell(OV.ng,1);
OV.gene.min = zeros(OV.ng,1);
OV.gene.max = zeros(OV.ng,1);
report_problems = false;
for g=1:OV.ng
  chr1=[];start1=[];end1=[];
  tmp = regexp(OV.gene.desc{g},'(chr.*):(\d*)-(\d*) ','tokens');
  if ~isempty(tmp)
    chr1 = tmp{1}{1};
    start1 = str2double(tmp{1}{2});
    end1 = str2double(tmp{1}{3});
  end
  chr2=[];start2=[];end2=[];
  idx = find(OV.targ.gene==g);
  if ~isempty(idx)
    chr2 = unique(OV.targ.chr(idx));
    start2 = min(OV.targ.start(idx));
    end2 = max(OV.targ.end(idx));
  end
  if report_problems
    if length(chr2)>1
      fprintf('%s has more than one chromosome in OV\n', OV.gene.name{g});
      chr2
    end
    if isempty(chr1) && isempty(chr2)
%      fprintf('%s has no location information\n', OV.gene.name{g});
    end
    if ~isempty(chr1) && ~isempty(chr2)
      if ~strcmp(chr1,chr2)
        fprintf('%s has chromosome conflict between OV (%s) and Expr (%s)\n', OV.gene.name{g},chr2,chr1);
      end  
      if ~(start1<end2 && start2<end1)
        fprintf('%s has failed overlap between OV and Expr\n', OV.gene.name{g});
        start1,end1
        start2,end2
      end
    end
  end
  if ~isempty(chr1)
    OV.gene.chr{g} = chr1;
    OV.gene.min(g) = start1;
    OV.gene.max(g) = end1;
    idx = find(strcmp(Cyto.chr,chr1)&Cyto.end>start1,1);
    if ~isempty(idx), OV.gene.band{g} = Cyto.band{idx}; end
  elseif ~isempty(chr2)
    if iscell(chr2), chr2=chr2{1}; end
    OV.gene.chr{g} = chr2;
    OV.gene.min(g) = start2;
    OV.gene.max(g) = end2;
    idx = find(strcmp(Cyto.chr,chr2)&Cyto.end>start2,1);
    if ~isempty(idx), OV.gene.band{g} = Cyto.band{idx}; end
  end   
end

% load Entrez GeneIDs and long names from Hugo database

H = load_struct('/xchip/tcga/gbm/analysis/lawrence/db/hugo.txt');
A = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');
H.gene = apply_aliases(H.ApprovedSymbol,A);
idx = listmap(OV.gene.name,H.gene);
i1 = find(~isnan(idx));
i2 = idx(i1);
OV.gene.entrez = cell(OV.ng,1);
OV.gene.entrez(i1) = H.EntrezGeneID(i2);
OV.gene.longname = cell(OV.ng,1);
OV.gene.longname(i1) = H.ApprovedName(i2);


%%%
%%%  make note of which genes were sequenced on which platforms
%%%

c6klist = '/xchip/tcga_scratch/lawrence/capture/tcga_6k_genes.targets.interval_list_NUMSONLY_WITHGENES.txt';
c2klist = '/xchip/tcga_scratch/lawrence/capture/cancer_2000gene_shift170.targets.interval_list_NUMSONLY_WITHGENES.txt';
c6k = load_struct(c6klist,'%s%s%s%s',0); c6k_genes = unique(c6k.col1);
c2k = load_struct(c2klist,'%s%s%s%s',0); c2k_genes = unique(c2k.col1);
idx = listmap(c6k_genes, OV.gene.name);
OV.gene.included_in_c6k = false(slength(OV.gene),1);
OV.gene.included_in_c6k(idx(~isnan(idx))) = true;
idx = listmap(c2k_genes, OV.gene.name);
OV.gene.included_in_c2k = false(slength(OV.gene),1);
OV.gene.included_in_c2k(idx(~isnan(idx))) = true;


%%%
%%%  load cluster assignments if available
%%%

if ~isempty(P.cluster_assignments_file)
  tmp = load_struct(P.cluster_assignments_file,'%s%s',0);
  OV.patient.cluster = map_across(OV.patient.short,regexprep(tmp.col1,'^TCGA-..-(....)$','$1'),tmp.col2);
  OV.patient.cluster = str2double(OV.patient.cluster);
end

