function find_frequently_mutated_genesets(MM,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'geneset_collection_file','/xchip/tcga/gbm/analysis/lawrence/db/gsea_canonical_pathway_plus_HMT67.txt');
P=impose_default_value(P,'freqgenesets_report_filename','*required*');

fprintf('Analyzing frequently mutated genesets\n');

M = MM.mut;
if isnumeric(M.gene)
  M.gene = nansub(MM.gene.name,M.gene);
end
if isnumeric(M.patient)
  M.patient = nansub(MM.patient.name,M.patient);
end
M = reorder_struct(M,grep('Missense|Nonsense|Splice_site',M.type,1));
g = unique(M.gene);
PW = which_pathways_include_these_genes(g,P);
npw = slength(PW);
PW.overlap_muts = repmat({''},npw,1);
PW.num_genes = zeros(npw,1);
PW.num_samples = zeros(npw,1);
PW.num_muts = zeros(npw,1);
for i=1:npw
  s = {}; m = [];
  for g=1:length(PW.overlap_genes{i})
    idx = find(strcmp(PW.overlap_genes{i}{g},M.gene));
    for j=1:length(idx), k=idx(j);
      m = [m; k]; s = [s; M.patient(k)];
      PW.overlap_muts{i} = [PW.overlap_muts{i} M.gene{k} '(' M.type{k} ':' M.patient{k} ') '];
    end
  end
  PW.num_genes(i) = length(PW.genes{i});
  PW.num_muts(i) = length(unique(m));
  PW.num_samples(i) = length(unique(s));
end
PW.overlap_muts = fillblanks(PW.overlap_muts,'(none)');
PW2 = PW;
for i=1:npw, PW2.genes{i} = concat(PW.genes{i},' '); end
PW2 = rmfield(PW2,{'num_overlap','overlap_genes'});
PW2 = sort_struct(PW2,{'num_samples','num_muts','num_genes'},[-1 -1 1]);
fprintf('Outputting %s\n',P.freqgenesets_report_filename);
save_struct(PW2,P.freqgenesets_report_filename);
