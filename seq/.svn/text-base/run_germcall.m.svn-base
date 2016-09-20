function run_germcall(sample,tn)

if ~exist('tn','var'), tn='n'; end

if tn(1)=='t', tn='tumor';
elseif tn(1)=='n', tn='normal';
else error('tn should be t or n');
end

boomdir = ['/xchip/tcga_scratch/ng/' sample '/wgs/boom/' tn '.boom'];
outdir = ['/xchip/tcga_scratch/ng/' sample '/wgs/' tn '_germcall'];

scatter_germcall(boomdir,outdir,sample);

