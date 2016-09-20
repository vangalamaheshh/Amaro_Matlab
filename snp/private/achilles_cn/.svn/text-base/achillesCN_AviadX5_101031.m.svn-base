% variation on script for comparison of Achilles and SCNA data
%
%   20010-10-31 give Aviad's data 10 more test runs for CDKN2A and
%   amplifications, 3 more for deletions to investigate the variance
%

%% input identifier for this iteration of the analysis

DATA_ID = 'Aviad100623';

input_dir = '/xchip/gistic/Achilles_correlates/achilles_data/';
output_dir = '/xchip/gistic/Achilles_correlates/Aviad_101028/';

% Aviad's gene score data file (note: not PMAD)
gscore_file = [input_dir 'Aviad/ach.v100623.Gs.txt'];
clmap_file = [input_dir 'Aviad/CLE-Aviad_CLmatch_100623.tab'];
mismatch_file = [output_dir DATA_ID '.unmatched_cell_lines.tab'];

% path to reference genome matlab file
ref_gene_file = ['/xchip/gistic/variables/hg18/' ...
                 'hg18_with_miR_20091116.mat'];

%% run parameters
validation_genes = true; % set to add fake genes with known outcomes to input

%% parameters

remove_X = true;
remove_cnvs = true;

% thresholds for calling amplifications and deletions
thresh = struct;
%...any
amp_thresh = log2(3)-1;
del_thresh = log2(1.5)-1;
%...none (control)
thresh.noamp = log2(2.3)-1;
thresh.nodel = log2(1.7)-1;
%...broad
thresh.broamp = log2(3)-1;
thresh.brodel = log2(1.7)-1;
%...focal
thresh.focamp = log2(3)-1;
thresh.focdel = log2(1.7)-1;

% thresholds for calling broad versus focal events
thresh.broad_fraction = 0.95;
thresh.focal_fraction = 0.95;

options.upper_amp_cap = -.5;
options.lower_del_cap = .5;

% options for permutation statistics
permtest_options = struct;
permtest_options.min_group = 7;          % minimal number of samples for each group
permtest_options.permType = 'binomial';  % permutation type

%% CELL LINE ENCYCLOPEDIA PROCESSING PARAMETERS

% CLE input files
cle_data_dir = '/xchip/gistic/CLE_projects/CCLE_data/';
% segmented copy number data
cle_input.segfile = [cle_data_dir 'CCLE_copynumber_2010-09-29.seg.txt'];
% cell line sample information
cle_input.sampinfo_file = [cle_data_dir 'CCLE_sample_info_file_2010-09-29.txt'];

% genomic locations of copy-number probes
cle_input.markersfile = ['/xchip/tcga_scratch/gsaksena/CancerGenomeAnalysisData/' ...
               'trunk/markerfiles/gistic_ovarian/broad.probes.txt'];
% regions that frequently have germ-line variations
cle_input.cnv_file = ['/xchip/tcga/gbm/analysis/mokelly/' ...
            '080429_convert_CNV_to_BED/' ...
            'CNV.verified_080606.combined.Ovarian.081202.txt'];

% CLE cache files
cle_cache.Qs_file = [cle_data_dir 'obsolete!'];
cle_cache.D_file = [cle_data_dir 'CLE_D.all.101012.smooth_10.mat'];
cle_cache.G_file = [cle_data_dir 'CLE_G.minmax.101012.symb.mat'];

% CLE processing options
cle_options = struct;
cle_options.use_segarray = true;
cle_options.join_segment_size = 10;
cle_options.remove_X = true;
cle_options.remove_nans = true;
cle_options.broad_length_cutoff = thresh.broad_fraction;
cle_options.focal_length_cutoff = thresh.focal_fraction;
cle_options.cn_cap = 1.5;

cle_options.gene_minmax = true;

%%% BEGIN PROCESSING %%%

%clear all
dbstop if error
set_verbose_level(30);

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
    
%% read hairpin data
verbose('=> Reading gene score data from ''%s''.',10,gscore_file);
AC = read_ataris_scores(gscore_file);

%% load reference genome into 'rg' and 'cyto'
load(ref_gene_file);
rg = rg([rg.chrn] < 100); % clean out haplotype variation entries in refgene
%% read or create CLE gene-level data

[G SI] = process_cle_data(cle_input,cle_cache,cle_options,rg,cyto);
%[G SI CL Qs] = process_cle_data(cle_input,cle_cache,cle_options,rg,cyto);

% match sample data (TODO: make into function)
[acx,gcx] = match_cell_lines(AC,G,SI,clmap_file,mismatch_file);
%!checkload_datasets;

% reduce CN data to matching samples
if exist('G','var')
    G = reorder_D_cols(G,gcx);
end
if exist('CL','var')
    CL = reorder_D_cols(CL,gcx);
end
if exist('Qs','var')
    Qs = subselect_Q(Qs,gcx);
end

% results structure
R = struct;

% match hairpins
[~,agx,cgx] = match_string_sets_hash(AC.geneID,G.geneID);
gmatch = sparse(agx,cgx,true,length(AC.geneID),length(G.geneID));
R.acgene = agx;
R.clgene = cgx;

% generate sample classes from cell-of-origin
cootype = {G.sis.siteprimary};
ucooty = unique(cootype);

classes = cell(1,length(ucooty));
for c=1:length(ucooty)
    classes{c} = strmatch(ucooty{c},cootype);
end

%% perform permutation statistics

% make amplification groups
anyamps = G.datmin >= amp_thresh;
noamps = G.datmax < thresh.noamp;
% define deletion groups
anydels = G.datmax <= del_thresh & G.datmin > log2(0.5)-1;
nodels = G.datmin > thresh.nodel;

scores = AC.dat(agx,acx);

N = 5; % number of replicates

xxx = rand(7,7); % roll the deterministic dice 

%% (1) CDKN2A dependence

CDKN2A_dependent = {...
    'COLO741_SKIN',...
    'COV504_OVARY',...
    'SW1990_PANCREAS',...
    'KYSE450_OESOPHAGUS',...
    'HEYA8_OVARY',...
    'LN229_CENTRAL_NERVOUS_SYSTEM',...
    'BXPC3_PANCREAS',...
    'SU8686_PANCREAS',...
    'KP4_PANCREAS',...
    'MM1S_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE',...
    'RMGI_OVARY' ...
};
[~,~,depmatch] = match_string_sets_hash(CDKN2A_dependent,G.sdesc);
dep = false(size(anyamps));
dep(:,depmatch) = true;

CDKN2A_independent = {...
    'CAOV4_OVARY',...
    'FUOV1_OVARY',...
    'HT29_LARGE_INTESTINE',...
    'HS766T_PANCREAS',...
    'OE33_OESOPHAGUS',...
    'NIHOVCAR3_OVARY',...
    'EFO21_OVARY',...
    'MKN7_STOMACH',...
    'OV90_OVARY',...
    'EFO27_OVARY',...
    'KYSE150_OESOPHAGUS',...
    'CFPAC1_PANCREAS',...
    'NCIH82_LUNG',...
    'GCIY_STOMACH',...
    'SJSA1_BONE',...
    'CAOV3_OVARY',...
    'OVCAR4_OVARY',...
    'HPAFII_PANCREAS' ...
};
[~,~,indmatch] = match_string_sets_hash(CDKN2A_independent,G.sdesc);
ind = false(size(noamps));
ind(:,indmatch) = true;

permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 7;
verbose('Permutation tests on CDKN2A dependence',20);
for k = 1:N
    tic
    R = hairpin_permtest(@stat_diffmean,scores,dep(cgx,:),ind(cgx,:),classes,permtest_options);
    toc
    if k == 1
        Ps = nan(size(R.pvalueL,1),N,2);
        qs = nan(size(R.pvalueL,1),N,2);
    end
    Ps(:,k,1) = R.pvalueL;
    Ps(:,k,2) = R.pvalueR;
    qs(:,k,1) = R.qvalueL;
    qs(:,k,2) = R.qvalueR;
end
% save intermediate results
n = write_multihp_permstats_101031([output_dir DATA_ID '.CDKN2A.diffmeanX' num2str(N) '.LR.bin.101031.tab'],Ps,qs,gmatch,AC);
verbose('wrote %d unfiltered amplification permutation statistics across hairpins',20,n);

%% (1a) enumerate all 17325 permutations

permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 7;
permtest_options.max_exact = 25000;
verbose('Exact permutation tests on CDKN2A dependence',20);
tic
R.xdep = hairpin_permtest(@stat_diffmean,scores,dep(cgx,:),ind(cgx,:),classes,permtest_options);
toc
% save intermediate results
n = write_hp_permstats([output_dir DATA_ID '.CDKN2A.diffmean.LR.bin.exact.101029.tab'],R.xdep,gmatch,AC,G,cyto,'dep','ind');
verbose('wrote %d unfiltered amplification permutation statistics across hairpins',20,n);


%% (2) deletions
permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 7;
verbose('Permutation tests on deletions',20);
for k = 1:N
    tic
    R = hairpin_permtest(@stat_diffmean,scores,anydels(cgx,:),nodels(cgx,:),classes,permtest_options);
    toc
    if k == 1
        Ps = nan(size(R.pvalueL,1),N,2);
        qs = nan(size(R.pvalueL,1),N,2);
    end
    Ps(:,k,1) = R.pvalueL;
    Ps(:,k,2) = R.pvalueR;
    qs(:,k,1) = R.qvalueL;
    qs(:,k,2) = R.qvalueR;
end
% save intermediate results
n = write_multihp_permstats_101031([output_dir DATA_ID '.alldels.HP.diffmeanX' num2str(N) '.LR.bin.101031.tab'],Ps,qs,gmatch,AC);
verbose('wrote %d unfiltered deletion permutation statistics across hairpins',20,n);

%% (3) amplifications
permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 7;
verbose('Permutation tests on amplifications',20);
for k = 1:N
    tic
    R = hairpin_permtest(@stat_diffmean,scores,anyamps(cgx,:),noamps(cgx,:),classes,permtest_options);
    toc
    if k == 1
        Ps = nan(size(R.pvalueL,1),N,2);
        qs = nan(size(R.pvalueL,1),N,2);
    end
    Ps(:,k,1) = R.pvalueL;
    Ps(:,k,2) = R.pvalueR;
    qs(:,k,1) = R.qvalueL;
    qs(:,k,2) = R.qvalueR;
end
% save intermediate results
n = write_multihp_permstats_101031([output_dir DATA_ID '.allamps.diffmeanX' num2str(N) '.LR.bin.101031.tab'],Ps,qs,gmatch,AC);
verbose('wrote %d unfiltered amplification permutation statistics across hairpins',20,n);

