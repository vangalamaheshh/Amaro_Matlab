function gset = load_genesets(P,M)

if ~exist('P','var'), P=[]; end

if ischar(P)
  tmp = P;
  P = [];
  P.geneset_collection_file = tmp;
end

P=impose_default_value(P,'geneset_collection_file','/cga/tcga-gsc/home/lawrence/xchip_tcga_gbm_analysis_lawrence/db/gsea_canonical_pathway.txt');
P=impose_default_value(P,'geneset_exclude','CANCER|MELANOMA|LEUKEMIA|GLIOMA|SCLEROSIS|DISEASE|DIABETES|CARCINOMA');
P=impose_default_value(P,'geneset_analysis_excludes_genes',{});

% load genesets

gset = load_gsea_collection(P.geneset_collection_file,'remove_slashes');
ngset = slength(gset);


for i=1:ngset
  gset.ngenes(i,1) = length(gset.genes{i});
end




if ~exist('M','var'), return; end

% which genes to exclude

masked_gene_names = M.gene.name;
excl = listmap(P.geneset_analysis_excludes_genes,M.gene.name);
masked_gene_names(excl) = repmat({'***excluded***'},length(excl),1);

% map genesets to project genes

fprintf('Mapping %d genesets from %s to list of project genes.\n', ngset, P.geneset_collection_file);
gset.genenos = cell(ngset,1);
for i=1:ngset
  if ~mod(i,100), fprintf('%d/%d ',i,ngset); end
  gset.genenos{i} = listmap(gset.genes{i}, masked_gene_names);
  gset.genenos{i} = unique(gset.genenos{i}(~isnan(gset.genenos{i})));
end
fprintf('\n');

% pick genesets to use

nonempty = find(~cellfun('isempty',gset.genenos));
excl = grep(P.geneset_exclude,gset.name,1);
use = setdiff(nonempty,excl);
gset = reorder_struct(gset,use);
ngset = slength(gset);



