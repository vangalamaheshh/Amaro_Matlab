function do_cov2_wig(samples)

if ~iscell(samples), samples = {samples}; end

extract_from_wig(samples,['/xchip/tcga_scratch/lawrence/capture/' ...
                    'whole_exome_agilent_designed_120.targets.interval_list.GENESGC_noheader.txt'],...
  'we16k_genes_category_coverage.txt','/xchip/tcga_scratch/lawrence/db/context',1,4);
