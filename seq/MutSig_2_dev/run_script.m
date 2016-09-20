%% Load DLBCL data 

P=[];
P.build = 'hg19';
P.targlist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.txt';
P.categdir = '/xchip/cga1/lawrence/db/hg19/context65';
P.categfile = '/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_hg19_june2010_terr_only.wig';
P.catfile = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/categs_CpGtransit_otherCGtransit_CGtransver_ATmut.txt';
P.genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
P.summed_cov_track = '/xchip/cga1/firehose_output/DLBCL_Analysis/Individual_Set/An_DLBCL_20110815/mutsig_1.5/An_DLBCL_20110815.summed_coverage.fwb';
isetname = 'An_DLBCL_20110815';
indir = '/xchip/cga1/firehose_output/DLBCL_Analysis/Individual_Set/An_DLBCL_20110815/mutsig_1.5';
P.patlist = [indir '/' isetname '.patients.txt'];
P.mutfile = [indir '/' isetname '.maf'];
P.ignore_categ_list_mismatch = true;   % (because coverage file includes indel+null category)
P.ignore_coverage = true;
M_dlbcl = load_all_mutation_data3(P,'dlbcl');

%% Load Lung Data

P=[];
P.build = 'hg19';
P.targlist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.txt';
P.categdir = '/xchip/cga1/lawrence/db/hg19/context65';
P.categfile = '/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_hg19_june2010_terr_only.wig';
P.catfile = '/xchip/cga1/firehose_output/Melanoma_all/Individual_Set/PR_Melanoma_CIP_Capture/mutsig/PR_Melanoma_CIP_Capture.mutcategs.txt';
P.genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
P.summed_cov_track = '/xchip/cga1/firehose_output/Melanoma_all/Individual_Set/PR_Melanoma_CIP_Capture/mutsig/PR_Melanoma_CIP_Capture.summed_coverage.fwb';
isetname = 'PR_Melanoma_CIP_Capture';
indir = '/xchip/cga1/firehose_output/Melanoma_all/Individual_Set/PR_Melanoma_CIP_Capture/mutsig/';
P.patlist = [indir '/' isetname '.patients.txt'];
P.mutfile = [indir '/' isetname '.maf'];
P.ignore_categ_list_mismatch = true;   % (because coverage file includes indel+null category)
P.ignore_coverage = true;
M_mel = load_all_mutation_data3(P,'melanoma');
clear P 
P.out_file = '/xchip/cga1/petar/Mutsig_2_dev/Mutsig_Addon/final_melanoma_results.txt';
P.image_dir = '/xchip/cga1/petar/Mutsig_2_dev/Mutsig_Addon/images_melanoma';

%% Load CLL Data 
P=[];
P.build = 'hg18';
P.targlist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg18.txt';
P.categdir = '/xchip/cga1/lawrence/db/context65';
P.categfile = '/xchip/cga1/lawrence/db/context65/RefSeq_exons_hg18_june2010_terr_only.wig';
P.catfile = '/xchip/cga1/firehose_output/CLL_paper/Individual_Set/An_CLL_88WE_3WGS_II/mutsig/An_CLL_88WE_3WGS_II.mutcategs.txt';
P.genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
P.summed_cov_track = '/xchip/cga1/firehose_output/CLL_paper/Individual_Set/An_CLL_88WE_3WGS_II/mutsig/An_CLL_88WE_3WGS_II.summed_coverage.fwb';
isetname = 'An_CLL_88WE_3WGS_II';
indir = '/xchip/cga1/firehose_output/CLL_paper/Individual_Set/An_CLL_88WE_3WGS_II/mutsig/';
P.patlist = [indir '/' isetname '.patients.txt'];
P.mutfile = ['/xchip/cga1/petar/CLL/An_CLL_88WE_3WGS_II.final_filtered_v4.maf'];
P.ignore_categ_list_mismatch = true;   % (because coverage file includes indel+null category)
P.ignore_coverage = true;
M_cll = load_all_mutation_data3(P,'cll');

M = {M_cll, M_dlbcl};

call_mutsig_cluster(M);