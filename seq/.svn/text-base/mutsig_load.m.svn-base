function M = mutsig_load(P,isetname)

if nargin<1, error('need P'); end
if ~exist('isetname','var'), isetname = 'iset'; end

if ischar(P)
  M = load_all_mutation_data(P);
  P=[];
else
  P = impose_default_value(P,'isetname',isetname);
  P = impose_default_value(P,'build','hg19');
  P = impose_default_value(P,'mutfile','*required*');
  P = impose_default_value(P,'covfile','*required*');
  P = impose_default_value(P,'BMR_covariates',  {...
      '/xchip/cga1/lawrence/db/hg19/covariates/RT_extra_median.txt';...
      '/xchip/cga1/lawrence/db/hg19/covariates/GC_300Kb.txt';...
      '/xchip/cga1/lawrence/db/hg19/covariates/gene_density_5Mb.txt';...
      '/xchip/cga1/lawrence/db/hg19/covariates/compartment.txt';...
      '/xchip/cga1/lawrence/db/hg19/covariates/log10_RNA_Seq_LUSC.txt';...
  });
  P = impose_default_value(P,'exclude_genes_equal_to_zero', ...
      '/xchip/cga1/lawrence/db/hg19/covariates/RNA_Seq_LUSC.txt');
  M = load_all_mutation_data(P);
end

P = impose_default_value(P,'precomputed_bagels', ...
  '/xchip/cga1/lawrence/db/hg19/covariates/bagels.v1.mat');

if ~isempty(P.precomputed_bagels) && ~strcmpi(P.precomputed_bagels,'none')
  fprintf('Loading precomputed bagels.\n');
  load(P.precomputed_bagels,'bagels','gname');
  map1 = listmap(M.gene.name,gname);
  map2 = listmap(gname,M.gene.name);
  bagels2 = bagels;
  for i=1:size(bagels2,2), bagels2(:,i) = nansub(map2,bagels(:,i)); end
  M.gene.bagel = nansub(bagels2,map1);
end


