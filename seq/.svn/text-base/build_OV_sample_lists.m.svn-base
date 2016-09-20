function [C2K C6K WE WGS] = build_OV_sample_lists(purpose)

if ~exist('purpose','var')
  purpose = 'mutation_calling_coverage';
end

C2K = []; C2K.sample.dir = {}; C6K = C2K; WE = C2K; WGS = C2K;

C2K.captype = 'C2K'; C6K.captype = 'C6K'; WE.captype = 'WE'; WGS.captype = 'WE'; %(sic)
C2K.file.targ = ['/xchip/tcga_scratch/lawrence/capture/' ...
           'cancer_2000gene_shift170.targets.interval_list.GENESGC.txt'];
C6K.file.targ = ['/xchip/tcga_scratch/lawrence/capture/' ...
   'tcga_6k_genes.targets.interval_list.GENESGC.txt'];
WE.file.targ = ['/xchip/tcga_scratch/lawrence/capture/' ...
   'whole_exome_refseq_coding.targets.interval_list.GENESGC.txt'];
WGS.file.targ = WE.file.targ;

categdir ='/xchip/tcga_scratch/lawrence/db/context'; numcategs = 4;
d = dir('/xchip/tcga_scratch/lawrence/ov/*');
idx2 = 1; idx6 = 1; idxwe = 1; idxwgs = 1;
for i=1:length(d)
  stem = ['/xchip/tcga_scratch/lawrence/ov/' d(i).name];
  t = dir([stem '/c2k/tumor_cbb/*.cbb']); n = dir([stem '/c2k/normal_cbb/*.cbb']);
  if length(t)>=24 & length(n)>=24
    C2K.sample.dir{idx2,1} = ['ov/' d(i).name '/c2k'];
    C2K.sample.tbam{idx2,1} = [stem '/c2k/tumor.bam'];
    C2K.sample.nbam{idx2,1} = [stem '/c2k/normal.bam'];
    idx2=idx2+1;
  end
  t = dir([stem '/c6k/tumor_cbb/*.cbb']); n = dir([stem '/c6k/normal_cbb/*.cbb']);
  if length(t)>=24 & length(n)>=24
    C6K.sample.dir{idx6,1} = ['ov/' d(i).name '/c6k'];
    C6K.sample.tbam{idx6,1} = [stem '/c6k/tumor.bam'];
    C6K.sample.nbam{idx6,1} = [stem '/c6k/normal.bam'];
    idx6=idx6+1;
  end
  t = dir([stem '/we/tumor_cbb/*.cbb']); n = dir([stem '/we/normal_cbb/*.cbb']);
  if length(t)>=24 & length(n)>=24
    WE.sample.dir{idxwe,1} = ['ov/' d(i).name '/we'];
    WE.sample.tbam{idxwe,1} = [stem '/we/tumor.bam'];
    WE.sample.nbam{idxwe,1} = [stem '/we/normal.bam'];
    idxwe=idxwe+1;
  end
  t = dir([stem '/wgs/tumor_cbb/*.cbb']); n = dir([stem '/wgs/normal_cbb/*.cbb']);
  if length(t)>=24 & length(n)>=24
    WGS.sample.dir{idxwgs,1} = ['ov/' d(i).name '/wgs'];
    WGS.sample.tbam{idxwgs,1} = [stem '/wgs/tumor.bam'];
    WGS.sample.nbam{idxwgs,1} = [stem '/wgs/normal.bam'];
    idxwgs=idxwgs+1;
  end
end

C2K.sample.short = regexprep(C2K.sample.dir,'ov/(....)/.{2,3}','$1');
C6K.sample.short = regexprep(C6K.sample.dir,'ov/(....)/.{2,3}','$1');
WE.sample.short = regexprep(WE.sample.dir,'ov/(....)/.{2,3}','$1');
WGS.sample.short = regexprep(WGS.sample.dir,'ov/(....)/.{2,3}','$1');

C2K.sample.med = regexprep(C2K.sample.short,'(.*)','OV-$1');
C6K.sample.med = regexprep(C6K.sample.short,'(.*)','OV-$1');
WE.sample.med = regexprep(WE.sample.short,'(.*)','OV-$1');
WGS.sample.med = regexprep(WGS.sample.short,'(.*)','OV-$1');


C2K.bad = load_struct('/xchip/tcga_scratch/lawrence/ov/mixedup_samples.txt','%s',0);
C2K.bad = unique(regexprep(C2K.bad.col1,'TCGA-..-(....)-...-...','$1'));
C2K.bad = [C2K.bad; '0931'];  % 0931 is fallopian
C2K.bad = [C2K.bad; '0757'];  % 0757 has almost zero normal coverage
C2K.sample = reorder_struct(C2K.sample,~ismember(C2K.sample.short,C2K.bad));

C6K.bad = {'1356';'1357'};  % bams for 1356 and 1357 have disappeared
C6K.sample = reorder_struct(C6K.sample,~ismember(C6K.sample.short,C6K.bad));

WE.bad = {'0188'};  % 0188 is a GBM
WE.bad = [ WE.bad; '0924'] ; % 0924 has almost zero on-target tumor coverage
WE.sample = reorder_struct(WE.sample,~ismember(WE.sample.short,WE.bad));

if strcmp(purpose,'mutation_calling_coverage')
  C2K.file.categdir = categdir; C2K.ns = slength(C2K.sample); C2K.ncat = numcategs;
  C6K.file.categdir = categdir; C6K.ns = slength(C6K.sample); C6K.ncat = numcategs;
  WE.file.categdir = categdir; WE.ns = slength(WE.sample); WE.ncat = numcategs;
  WGS.file.categdir = categdir; WGS.ns = slength(WGS.sample); WGS.ncat = numcategs;

  C2K.file.cov = 'c2k_genes_category_coverage.txt';
  C6K.file.cov = 'c6k_genes_category_coverage.txt';
  WE.file.cov = 'we_genes_category_coverage.txt';
  WGS.file.cov = WE.file.cov;

elseif strcmp(purpose,'TN_coverage')
  C2K.file.categdir = []; C2K.ns = slength(C2K.sample); C2K.ncat = 1;
  C6K.file.categdir = []; C6K.ns = slength(C6K.sample); C6K.ncat = 1;
  WE.file.categdir = []; WE.ns = slength(WE.sample); WE.ncat = 1;
  WGS.file.categdir = []; WGS.ns = slength(WGS.sample); WGS.ncat = 1;

  C2K.file.cov = 'c2k_TN_coverage.txt';
  C6K.file.cov = 'c6k_TN_coverage.txt';
  WE.file.cov = 'we_TN_coverage.txt';
  WGS.file.cov = WE.file.cov;

else
  error('Unknown purpose');
end


