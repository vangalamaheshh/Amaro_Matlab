function do_cov3(samples)

if ~iscell(samples), samples = {samples}; end

reglist = '/xchip/tcga_scratch/lawrence/db/genome.txt';
outfile = 'genome_allcateg_coverage.txt';
categdir = '/xchip/tcga_scratch/lawrence/db/allcateg';
ncategs = 1040;
extract_from_cbb(samples,reglist,outfile,categdir,ncategs);
