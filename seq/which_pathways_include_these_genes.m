function PW = which_pathways_include_these_genes(g,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'geneset_collection_file','/xchip/tcga/gbm/analysis/lawrence/db/gsea_canonical_pathway.txt');
% alternative: '/xchip/tcga/gbm/analysis/lawrence/db/gsea_kegg.txt'
P=impose_default_value(P,'geneset_exclude','CANCER|MELANOMA|LEUKEMIA|GLIOMA|SCLEROSIS|DISEASE|DIABETES|CARCINOMA');
P=impose_default_value(P,'geneset_analysis_excludes_genes',{});

% load genesets

gset = load_gsea_collection(P.geneset_collection_file,'remove_slashes');
ngset = slength(gset);

% find pathways

PW = gset;
PW.num_overlap = zeros(slength(PW),1);
PW.overlap_genes = cell(slength(PW),1);

for i=1:slength(PW)
  x = intersect(PW.genes{i},g);
  PW.num_overlap(i) = length(x);
  PW.overlap_genes{i} = x;
end
