function do_cov2(samples)

if ~iscell(samples), samples = {samples}; end

extract_from_cbb(samples,'/xchip/tcga_scratch/lawrence/capture/whole_exome_refseq_coding.targets.interval_list.GENESGC.txt',...
  'we_genes_category_coverage.txt','/xchip/tcga_scratch/lawrence/db/context',4);
