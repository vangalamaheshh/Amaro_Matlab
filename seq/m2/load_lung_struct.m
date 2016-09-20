function [ M ] = load_lung_struct()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
P=[];
P.build = 'hg19';
P.targlist = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.txt';
P.categdir = '/xchip/cga1/lawrence/db/hg19/context65';
P.categfile = '/xchip/cga1/lawrence/db/hg19/context65/RefSeq_exons_hg19_june2010_terr_only.wig';
P.catfile = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/categs_CpGtransit_otherCGtransit_CGtransver_ATmut_Null.txt';
P.genelist = '/xchip/cga1/lawrence/capture/RefSeq_exons_hg19_june2010.genelist.txt';

isetname = 'An_TCGA_LUSC_Capture_101118_noHM';
P.patlist = '/xchip/cga/gdac-prod/genepattern/jobResults/598252/An_TCGA_LUSC_Capture_101118_NoHM.patients.txt';
P.covfile = '/xchip/cga/gdac-prod/genepattern/jobResults/598253/An_TCGA_LUSC_Capture_101118_NoHM.coverage.mat';
P.mutfile = '/xchip/cga/gdac-prod/genepattern/jobResults/598252/An_TCGA_LUSC_Capture_101118_NoHM.maf';

M = load_all_mutation_data2(P);


end

