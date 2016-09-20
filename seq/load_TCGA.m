function TCGA = load_TCGA(P)
%
% load_TCGA(P)
%
% P can be a parameter structure
% or the name of a dataset:
%   full, phase1_72, phase1_84, phase1_91, latest_154, latest_139, ova
%

default_dataset = 'nov';

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
case 'nov'
  P.patient_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/154patients.txt';
  P.gene_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/1325genes.txt';
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20081113/tcga_20081113.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20081113/tcga_20081113.cov';
case 'full'
  P.mut_file='/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080812/tcga_20080819b.mut';
  P.cov_file='/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080812/TCGA_context_20080813.txt';
  P.gene_file='/xchip/tcga/gbm/analysis/lawrence/tcga/allgenes.txt';
  P.patient_file='/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/allpatients.txt';
case 'phase1_91'
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080702/tcga_20080711.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080502/context_20080606.txt';
  P.gene_file='/xchip/tcga/gbm/analysis/lawrence/tcga/601genes.txt';
  P.patient_file='/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/91patients.txt';
case 'phase1_84'
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080702/tcga_20080711.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080502/context_20080606.txt';
  P.gene_file='/xchip/tcga/gbm/analysis/lawrence/tcga/601genes.txt';
  P.patient_file='/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/84patients.txt';
case 'phase1_72'
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080702/tcga_20080711.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080502/context_20080606.txt';
  P.gene_file='/xchip/tcga/gbm/analysis/lawrence/tcga/601genes.txt';
  P.patient_file='/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/72patients.txt';
case 'oct_154'
  P.patient_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/154patients.txt';
  P.gene_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/1325genes.txt';
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20081016/tcga_20081016c.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080812/TCGA_context_20080813.txt';
case 'oct_139'
  P.patient_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/139patients.txt';
  P.gene_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/1325genes.txt';
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20081016/tcga_20081016c.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20080812/TCGA_context_20080813.txt';
case 'ova'
  P.patient_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081029/24patients.txt';
  P.gene_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/1325genes.txt';
  P.mut_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081126/maf/tcga_ova_20081202b.mut';
  P.cov_file = '/xchip/tcga/gbm/analysis/lawrence/tcga/ova/20081126/cov/tcga_ova_20081202.cov';
otherwise
  error('Unknown dataset');
end, end



% LOAD MUTATION DATA

TCGA = load_mutdata3(P);

P=impose_default_value(P,'use_validated_or_verified_mutations_only',true);
P=impose_default_value(P,'but_use_all_silent_mutations', true);
TCGA = choose_mutations(TCGA,P);
TCGA = build_n_and_N_tables(TCGA);
TCGA = process_3N_breakdown_stats(TCGA);

TCGA.gene.desc = repmat({''},TCGA.ng,1);
TCGA.patient.have_seq = squeeze(sum(sum(TCGA.N_cov(:,TCGA.TOT,:),1),2))>0;
TCGA.gene.have_seq = squeeze(sum(sum(TCGA.N_cov(:,TCGA.TOT,:),3),2))>0;

% LOAD ALIASES

aliases = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');

if isempty(grep('ova',TCGA.file.mut))

% LOAD EXPRESSION DATA

% E=read_R_expression_file('/xchip/tcga/gbm/analysis/roel/BLU/Originals/mainpaper/TCGA_206_median.txt');
%      copied to my directory as expr.txt
%
% > cat expr.txt | sed 's/\bNA\b/NaN/g' > expr_modified.txt

fprintf('Loading expression data\n');
TCGA.file.expr = '/xchip/tcga/gbm/analysis/roel/BLU/Originals/mainpaper/TCGA_206_median.txt';
ddir = '/xchip/tcga/gbm/analysis/lawrence/genefig/';
tmp = read_table([ddir 'expr_modified.txt'],['%s' repmat('%f',1,206)],char(9),1,'whitespace','\b\r');
tmp.gene.name = apply_aliases(tmp.dat{1},aliases);
tmp.patient.code = tmp.headers{1}(2:end)';
tmp.patient.name = regexprep(tmp.patient.code,'-...-..$','');
tmp.data = cat(2,tmp.dat{2:end});

% reorder to match TCGA patient and gene lists
gidx = listmap(TCGA.gene.name,tmp.gene.name);
pidx = listmap(TCGA.patient.name,tmp.patient.name);
tmp.data(:,end+1)=NaN;
tmp.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(tmp.data,1);
pidx(isnan(pidx))=size(tmp.data,2);
TCGA.expr = tmp.data(gidx,pidx);
TCGA.patient.have_expr = ~all(isnan(TCGA.expr),1)';
TCGA.gene.have_expr = ~all(isnan(TCGA.expr),2);

% LOAD COPY-NUMBER DATA

% C=read_eisen_dat('/xchip/tcga/gbm/analysis/mokelly/080512_call_genes/summarized_genes_spreadsheet.080612.txt');
%      copied to my directory as cn.txt

% full gene list in summarized_genes_spreadsheet.080630.txt
% C=read_eisen_dat('/xchip/tcga/gbm/analysis/mokelly/080512_call_genes/summarized_genes_spreadsheet.080630.txt');
%      copied to my directory as cn_full.txt

% > cat cn_full.txt | sed 's/\bInconsistent\b/Inf/g' > cn_full_modified.txt

fprintf('Loading copy number data\n');
TCGA.file.cn = '/xchip/tcga/gbm/analysis/mokelly/080512_call_genes/summarized_genes_spreadsheet.080630.txt';
ddir = '/xchip/tcga/gbm/analysis/lawrence/genefig/';
tmp = read_table([ddir 'cn_full_modified.txt'],['%s%s' repmat('%f',1,206)],char(9),1,'whitespace','\b\r');
tmp.gene.name = apply_aliases(tmp.dat{1},aliases);
tmp.gene.desc = tmp.dat{2};
tmp.patient.name = tmp.headers{1}(3:end)';
tmp.data = cat(2,tmp.dat{3:end});
tmp.data(isinf(tmp.data))=NaN;   % Inf->Nan

% reorder to match TCGA patient and gene lists
gidx = listmap(TCGA.gene.name,tmp.gene.name);
pidx = listmap(TCGA.patient.name,tmp.patient.name);
tmp.data(:,end+1)=NaN;
tmp.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(tmp.data,1);
pidx(isnan(pidx))=size(tmp.data,2);
TCGA.cn = tmp.data(gidx,pidx);
TCGA.patient.have_cn = ~all(isnan(TCGA.cn),1)';
TCGA.gene.have_cn = ~all(isnan(TCGA.cn),2);

tmp.gene.desc = [tmp.gene.desc; {''}];
TCGA.gene.desc = tmp.gene.desc(gidx);

% LOAD METHYLATION DATA

% -->orig file is at
%       /xchip/tcga/gbm/visualization/tcga_20080815/sorted.methylation.affy.080819.igv.txt
% --> copied it to my own directory as "methyl.txt"
% --> removed first line
% > tail -1506 methyl.txt > methyl_modified.txt
% --> changed N/A to NaN
% > cat methyl_modified.txt | sed 's/\bN\/A\b/NaN/g' > methyl_modified2.txt

TCGA.file.methyl = '/xchip/tcga/gbm/visualization/tcga_20080815/sorted.methylation.affy.080819.igv.txt';
ddir = '/xchip/tcga/gbm/analysis/lawrence/genefig/';
tmp = read_table([ddir 'methyl_modified2.txt'],['%s%f%f%s' repmat('%f',1,265)],char(9),1,'whitespace','\b\r');
tmp.gene.name = apply_aliases(regexprep(tmp.dat{4},'_.*',''),aliases);
tmp.patient.code = tmp.headers{1}(5:end)';
tmp.patient.name = regexprep(tmp.patient.code,'-...-...-....-..$','');
tmp.data = cat(2,tmp.dat{5:end});

% for each gene, take average methylation across probes for that gene
[g gi gj] = unique(tmp.gene.name);
keep = true(size(tmp.data,1),1);
for i=1:length(g)
  idx = find(gj==i);
  tmp.data(idx(1),:) = mean(tmp.data(idx,:));
  keep(idx(2:end)) = false;
end
tmp.data = tmp.data(keep,:);
tmp.gene = reorder_struct(tmp.gene,keep);

% reorder to match TCGA patient and gene lists
gidx = listmap(TCGA.gene.name,tmp.gene.name);
pidx = listmap(TCGA.patient.name,tmp.patient.name);
tmp.data(:,end+1)=NaN;
tmp.data(end+1,:)=NaN;
gidx(isnan(gidx))=size(tmp.data,1);
pidx(isnan(pidx))=size(tmp.data,2);
TCGA.methyl = tmp.data(gidx,pidx);
TCGA.patient.have_methyl = ~all(isnan(TCGA.methyl),1)';
TCGA.gene.have_methyl = ~all(isnan(TCGA.methyl),2);

end % if ~strcmp(tag,'ova')

% PARSE GENE LOCATIONS

fprintf('Parsing gene locations\n');

Cyto = load_struct('/xchip/tcga/gbm/analysis/lawrence/genome/hg18/cytoBand.txt','%s%f%f%s%s');

TCGA.gene.chr = cell(TCGA.ng,1);
TCGA.gene.band = cell(TCGA.ng,1);
TCGA.gene.min = zeros(TCGA.ng,1);
TCGA.gene.max = zeros(TCGA.ng,1);
report_problems = false;
for g=1:TCGA.ng
  chr1=[];start1=[];end1=[];
  tmp = regexp(TCGA.gene.desc{g},'(chr.*):(\d*)-(\d*) ','tokens');
  if ~isempty(tmp)
    chr1 = tmp{1}{1};
    start1 = str2double(tmp{1}{2});
    end1 = str2double(tmp{1}{3});
  end
  chr2=[];start2=[];end2=[];
  idx = find(TCGA.targ.gene==g);
  if ~isempty(idx)
    chr2 = unique(TCGA.targ.chr(idx));
    start2 = min(TCGA.targ.start(idx));
    end2 = max(TCGA.targ.end(idx));
  end
  if report_problems
    if length(chr2)>1
      fprintf('%s has more than one chromosome in TCGA\n', TCGA.gene.name{g});
      chr2
    end
    if isempty(chr1) && isempty(chr2)
%      fprintf('%s has no location information\n', TCGA.gene.name{g});
    end
    if ~isempty(chr1) && ~isempty(chr2)
      if ~strcmp(chr1,chr2)
        fprintf('%s has chromosome conflict between TCGA (%s) and Expr (%s)\n', TCGA.gene.name{g},chr2,chr1);
      end  
      if ~(start1<end2 && start2<end1)
        fprintf('%s has failed overlap between TCGA and Expr\n', TCGA.gene.name{g});
        start1,end1
        start2,end2
      end
    end
  end
  if ~isempty(chr1)
    TCGA.gene.chr{g} = chr1;
    TCGA.gene.min(g) = start1;
    TCGA.gene.max(g) = end1;
    idx = find(strcmp(Cyto.chr,chr1)&Cyto.end>start1,1);
    if ~isempty(idx), TCGA.gene.band{g} = Cyto.band{idx}; end
  elseif ~isempty(chr2)
    if iscell(chr2), chr2=chr2{1}; end
    TCGA.gene.chr{g} = chr2;
    TCGA.gene.min(g) = start2;
    TCGA.gene.max(g) = end2;
    idx = find(strcmp(Cyto.chr,chr2)&Cyto.end>start2,1);
    if ~isempty(idx), TCGA.gene.band{g} = Cyto.band{idx}; end
  end   
end

% load Entrez GeneIDs and long names from Hugo database

H = load_struct('/xchip/tcga/gbm/analysis/lawrence/db/hugo.txt');
A = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');
H.gene = apply_aliases(H.ApprovedSymbol,A);
idx = listmap(TCGA.gene.name,H.gene);
i1 = find(~isnan(idx));
i2 = idx(i1);
TCGA.gene.entrez = cell(TCGA.ng,1);
TCGA.gene.entrez(i1) = H.EntrezGeneID(i2);
TCGA.gene.longname = cell(TCGA.ng,1);
TCGA.gene.longname(i1) = H.ApprovedName(i2);
