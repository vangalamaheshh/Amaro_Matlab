function [tlanes nlanes] = get_lanecounts(sample)

tlanes = nan; nlanes = nan;

try
basedir = ['/xchip/tcga_scratch/lawrence/' sample];
t = load_lines([basedir '/tumor.bam.lanelist']);
n = load_lines([basedir '/normal.bam.lanelist']);
tlanes = length(t); nlanes = length(n);

catch me
  fprintf('UNABLE TO RETRIEVE LANECOUNTS FOR %s\n',sample);
end
