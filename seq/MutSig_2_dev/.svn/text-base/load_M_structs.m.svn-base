


P=[];
P.build = 'hg19';
P.targlist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.txt';
P.categdir = '/xchip/cga1/lawrence/db/hg19/context65';
P.categfile = '/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_hg19_june2010_terr_only.wig';
P.catfile = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/categs_CpGtransit_otherCGtransit_CGtransver_ATmut.txt';
P.genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
P.summed_cov_track = '/xchip/cga1/lawrence/nb/analysis/20100914/nb76/nb76.summed.somatic_coverage.fwb'; 
isetname = 'nb76';
indir = '/xchip/cga1/lawrence/nb/analysis/20100914/nb76';
P.patlist = [indir '/' isetname '.patients.txt'];
P.mutfile = [indir '/' isetname '.maf'];
P.ignore_categ_list_mismatch = true;   % (because coverage file includes indel+null category)
P.ignore_coverage = true;
M_nb = load_all_mutation_data3(P,'neuroblastoma');


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



P=[];
P.build = 'hg19';
P.targlist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg19.txt';
P.categdir = '/xchip/cga1/lawrence/db/hg19/context65';
P.categfile = '/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_hg19_june2010_terr_only.wig';
P.catfile = '/xchip/cga1/lawrence/lung/analysis/20101210/An_TCGA_LUSC_Capture_101118_noHM.mutcategs.null.txt';
P.genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
P.patlist = '/xchip/cga/gdac-prod/genepattern/jobResults/598252/An_TCGA_LUSC_Capture_101118_NoHM.patients.txt';
P.mutfile = '/xchip/cga/gdac-prod/genepattern/jobResults/598252/An_TCGA_LUSC_Capture_101118_NoHM.maf';
P.summed_cov_track = '/xchip/cga1/lawrence/lung/analysis/20101210/An_TCGA_LUSC_Capture_101118_NoHM.summed.somatic_coverage.fwb';
P.ignore_coverage = true;
M_lung = load_all_mutation_data3(P,'lung');


P=[];
P.build = 'hg18';
P.mutfile = '/xchip/cga1/lawrence/ov/analysis/20101110/ov316.maf';
P.patlist = '/xchip/cga1/lawrence/ov/analysis/20100818_3centers/ov316.patients.txt';
P.catfile = ['/xchip/cga2/lawrence/cga/trunk/matlab/seq/' ...
           'categs_CpGtransit_otherCGtransit_CGtransver_ATmut_IndelAndNull.txt'];
P.genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
P.ignore_categ_list_mismatch = true; % (because changed to indel+null)
P.targlist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_hg18.txt';
P.summed_cov_track = '/xchip/cga1/lawrence/ov/analysis/20100818_3centers/ov316.fwb/ov316.summed.somatic_coverage.fwb';
P.ignore_coverage = true;
M_ovarian = load_all_mutation_data3(P,'ovarian');

M = {M_nb, M_lung, M_ovarian};

M = {M_nb, M_ovarian}; 
