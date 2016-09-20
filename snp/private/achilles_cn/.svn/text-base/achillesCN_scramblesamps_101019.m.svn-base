% variation on script for comparison of Achilles and SCNA data
%
% - 2010-10-12 negative controls with scrambled samples for difference of
% mean and S2N runs with 2010.10.12 data
%
%% input identifier for this iteration of the analysis

DATA_ID = '20101012';

input_dir = '/xchip/gistic/Achilles_correlates/achilles_data/';
output_dir = '/xchip/gistic/Achilles_correlates/mingroup7_101012/';
output_dir2 = '/xchip/gistic/Achilles_correlates/S2N_101015/';

% achilles hairpin score data file (note: not PMAD)
hpscore_file = [input_dir '20101012_achilles2_modifiedscore.rnai.gct'];

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
verbose('=> Reading hairpin data from ''%s''.',10,hpscore_file);
AC = read_hairpin_data(hpscore_file);

%% load reference genome into 'rg' and 'cyto'
load(ref_gene_file);
rg = rg([rg.chrn] < 100); % clean out haplotype variation entries in refgene
%% read or create CLE gene-level data

[G SI] = process_cle_data(cle_input,cle_cache,cle_options,rg,cyto);
%[G SI CL Qs] = process_cle_data(cle_input,cle_cache,cle_options,rg,cyto);

% match sample data (TODO: make into function)
checkload_datasets;

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

% Z-score the scores for each cell line
means = mean(AC.dat,1);
stdevs = std(AC.dat,1);
ZAC = AC;
ZAC.dat = (AC.dat-repmat(means,size(AC.dat,1),1)) ./ repmat(stdevs,size(AC.dat,1),1);


% 19-Oct-2010 additional negative controls with scrambled samples
% for difference of means (101012) and S2N (101015)


%% (1) negative control - amplification with scrambled genes
SAC = ZAC;
SAC.dat = ZAC.dat(:,randperm(size(ZAC.dat,2)));
permtest_options.upperMeanLimit = Inf;
permtest_options.lowerMeanLimit = -Inf;
permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 6;
tic
R.nega = perform_perms_parallel(R,@(x,y) mean(x,1)-mean(y,1),SAC.dat(:,acx),anyamps,noamps,classes,permtest_options);
toc
% save intermediate results
n = write_hp_permstats([output_dir DATA_ID '.allamps.HP.scramblesamp.LR.bin.101019.tab'],R.nega,gmatch,SAC,G,cyto,'any amp','no amp');
verbose('wrote %d unfiltered deletion permutation statistics across hairpins',20,n);

%% (2) negative control - deletion with scrambled genes
SAC = ZAC;
SAC.dat = ZAC.dat(:,randperm(size(ZAC.dat,2)));
permtest_options.upperMeanLimit = Inf;
permtest_options.lowerMeanLimit = -Inf;
permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 6;
tic
R.negdel = perform_perms_parallel(R,@(x,y) mean(x,1)-mean(y,1),SAC.dat(:,acx),anydels,nodels,classes,permtest_options);
toc
% save intermediate results
n = write_hp_permstats([output_dir DATA_ID '.alldels.HP.scramblesamp.LR.bin.101019.tab'],R.negdel,gmatch,SAC,G,cyto,'any del','no del');

%% (3) negative control - amplification with scrambled genes
SAC = ZAC;
SAC.dat = ZAC.dat(:,randperm(size(ZAC.dat,2)));
permtest_options.upperMeanLimit = Inf;
permtest_options.lowerMeanLimit = -Inf;
permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 7;
permtest_options.nprocesses = 6;
tic
R.nega2 = perform_perms_parallel(R,@(x,y) (mean(x,1)-mean(y,1))./(std(x,1)+std(y,1)),SAC.dat(:,acx),anyamps,noamps,classes,permtest_options);
toc
% save intermediate results
n = write_hp_permstats([output_dir2 DATA_ID '.allamps.HP.scramblesamps.S2N.LR.bin.101019.tab'],R.nega2,gmatch,SAC,G,cyto,'any amp','no amp');
verbose('wrote %d unfiltered deletion permutation statistics across hairpins',20,n);

%% (4) negative control - deletion with scrambled genes
SAC = ZAC;
SAC.dat = ZAC.dat(:,randperm(size(ZAC.dat,2)));
permtest_options.upperMeanLimit = Inf;
permtest_options.lowerMeanLimit = -Inf;
permtest_options.permType = 'binomial';
permtest_options.tail = 'leftright';
permtest_options.min_group = 7;
permtest_options.nprocesses = 6;
tic
R.negd2 = perform_perms_parallel(R,@(x,y) (mean(x,1)-mean(y,1))./(std(x,1)+std(y,1)),SAC.dat(:,acx),anydels,nodels,classes,permtest_options);
toc
% save intermediate results
n = write_hp_permstats([output_dir2 DATA_ID '.alldels.HP.scramblesamps.S2N.LR.bin.101019.tab'],R.negd2,gmatch,SAC,G,cyto,'any del','no del');
verbose('wrote %d unfiltered deletion permutation statistics across hairpins',20,n);

