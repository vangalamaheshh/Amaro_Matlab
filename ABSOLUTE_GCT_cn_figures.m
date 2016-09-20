ABS_SEG=load_table('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/ABS_Segs/AggregatedGCTABS.seg'); ABS_SEG=rmfield(ABS_SEG,{'header','headline'});
samp1 = 'TGC-Testes_DFCI_18-TP-NT-SM-4PDE4-SM-4PDE6';
samp2 = 'TGC-Testes_DFCI_18-TP-NT-SM-4PDE5-SM-4PDE6';
load('rb_colomap.mat')
ABS_Table=load_table('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/AggregatedGCTABSTable.txt')
seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));
seg1.xstart=xhg19(seg1.Chromosome,seg1.Startbp); seg1.xend=xhg19(seg1.Chromosome,seg1.Endbp);
seg2.xstart=xhg19(seg2.Chromosome,seg2.Startbp); seg2.xend=xhg19(seg2.Chromosome,seg2.Endbp);
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));


samp1 = 'TGC-Testes_DFCI_7-TP-NT-SM-4PDDF-SM-4PDDJ';
samp2 = 'TGC-Testes_DFCI_7-TP-NT-SM-4PDDG-SM-4PDDJ';
samp3 = 'TGC-Testes_DFCI_7-TP-NT-SM-4PDDH-SM-4PDDJ';
seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));
seg3=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp3));
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp3)));
absolute_allelic_cn_plot( seg3,GD)




samp1 = 'TGC-Testes_DFCI_31-TP-NT-SM-4PDEY-SM-4PDF3';
samp2 = 'TGC-Testes_DFCI_31-TP-NT-SM-4PDEZ-SM-4PDF3';
samp3 = 'TGC-Testes_DFCI_31-TP-NT-SM-4PDF1-SM-4PDF3';
samp4 = 'TGC-Testes_DFCI_31-TP-NT-SM-4PDF2-SM-4PDF3';

seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));
seg3=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp3));
seg4=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp4));

GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp3)));
absolute_allelic_cn_plot( seg3,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp4)));
absolute_allelic_cn_plot( seg4,GD)

samp1='TGC-Testes_DFCI_9-TP-NT-SM-4PDDM-SM-4PDDO';
samp2='TGC-Testes_DFCI_9-TP-NT-SM-4PDDN-SM-4PDDO';
seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));

GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)

samp1='TGC-Testes_DFCI_18-TP-NT-SM-4PDE4-SM-4PDE6';
samp2='TGC-Testes_DFCI_18-TP-NT-SM-4PDE5-SM-4PDE6';
seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));

GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)


samp1='TGC-Testes_DFCI_4-TP-NT-SM-4PDD8-SM-4PDDA';
samp2='TGC-Testes_DFCI_4-TP-NT-SM-4PDD9-SM-4PDDA';
seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));

GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)

%ITGCN segementation files

samp1 = 'TGC-Testes_DFCI_55-ITGCN';
samp2 = 'TGC-Testes_DFCI_55-TP-NT-SM-7CMME-SM-7CMMF';
seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));

GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)


samp1 = 'TGC-Testes_DFCI_61-ITGCN';
samp2 = 'TGC-Testes_DFCI_61-TP-NT-SM-7CMMQ-SM-7CMMR';

seg1=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp1));
seg2=reorder_struct(ABS_SEG,ismember(ABS_SEG.sample,samp2));

GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp1)));
absolute_allelic_cn_plot( seg1,GD)
GD=ABS_Table.Genome_doublings(find(ismember(ABS_Table.sample,samp2)));
absolute_allelic_cn_plot( seg2,GD)
