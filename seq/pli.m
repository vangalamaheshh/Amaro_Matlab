function [t n] = pli(individual,chr,pos)
% plot primary data for a putative somatic insertion
% individual = e.g. 'MM-0309'

win = 10000;
sw = 300;
sw2 = 50;

t = bam;t.chr = chr; t.st = pos-win; t.en = pos+win;n = t;
t = t.loadbam(['/xchip/cga1/firehose_output/Individual/' individual '/wgs/bam/tumor.bam']);
n = n.loadbam(['/xchip/cga1/firehose_output/Individual/' individual '/wgs/bam/normal.bam']);
t = t.calc(sw,sw2,true); n = n.calc(sw,sw2,true);

plot_insertion_data(individual,chr,pos,win,t,n);
