function Mg = map_to_genesets(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'geneset_collection_file','/xchip/cga1/annotation/db/genesets/gsea_canonical_pathway_plus_HMT67.txt');
%P=impose_default_value(P,'geneset_collection_file','/xchip/tcga/gbm/analysis/lawrence/db/gsea_canonical_pathway.txt');
% alternative: '/xchip/tcga/gbm/analysis/lawrence/db/gsea_kegg.txt'
P=impose_default_value(P,'geneset_exclude','CANCER|MELANOMA|LEUKEMIA|GLIOMA|SCLEROSIS|DISEASE|DIABETES|CARCINOMA');
P=impose_default_value(P,'geneset_analysis_excludes_genes',{});

% load genesets

gset = load_gsea_collection(P.geneset_collection_file,'remove_slashes');
ngset = slength(gset);

% which genes to exclude

masked_gene_names = M.gene.name;
if ~isempty(P.geneset_analysis_excludes_genes)
  if ~iscell(P.geneset_analysis_excludes_genes)
    error('P.geneset_analysis_excludes_genes should be a cell array of strings (gene names to exclude)');
  end
  excl = listmap(P.geneset_analysis_excludes_genes,M.gene.name);
  masked_gene_names(excl) = repmat({'***excluded***'},length(excl),1);
%  fprintf('Excluding the following genes:\n');
%  for i=1:length(P.geneset_analysis_excludes_genes)
      %keyboard
     % fprintf('%s --> %s\n',P.geneset_analysis_excludes_genes{i},...
     %     decell(nansub({'found','NOT FOUND'},1+isnan(excl(i)))));
 % end
end

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

% build new geneset-centric data structure

Mg = M;
Mg = rename_field(Mg,'gene','gene_before_geneset_collapse');
Mg = rename_field(Mg,'ng','ng_before_geneset_collapse');
Mg.gene = [];
Mg.ng = ngset;
fprintf('Total of %d genesets contain at least one project gene.\n', Mg.ng);

% collapse N and n tables
flds = grep('^(N|n)_.*',fieldnames(M));
Mg = rmfield(Mg,flds);
for f=1:length(flds)
  old = getfield(M,flds{f});
  oldsize = size(old);
  if oldsize(1)~=M.ng, continue; end
  newsize = oldsize;
  newsize(1) = Mg.ng;
  new = zeros(newsize);
  for i=1:Mg.ng
    genes = gset.genenos{i};
    new(i,:,:,:,:) = sum(old(genes,:,:,:,:),1);
  end
  Mg = setfield(Mg,flds{f},new);
end

% tally hits in each gene
gset.muttally = cell(Mg.ng,1);
gset.mutnos = cell(Mg.ng,1);
for i=1:Mg.ng
  muttally = [];
  mutnos = [];
  genes = gset.genenos{i};
  for g=1:length(genes)
    n = sum(M.n_nonsilent(genes(g),M.TOT,:),3);
    if n   % if mutation(s) in this gene
      muttally = [muttally sprintf('%s(%d), ', M.gene.name{genes(g)}, n)];
      if isfield(M,'use_nonsilent')
        idx = M.use_nonsilent;
      elseif isfield(M.mut,'is_coding') && isfield(M.mut,'is_silent')
        idx = find(M.mut.is_coding & ~M.mut.is_silent);
      elseif isfield(M.mut,'type')
        idx = find(~strcmpi('silent',M.mut.type) && ~strcmpi('synonymous',M.mut.type)&&...
                   ~strcmpi('intron',M.mut.type) && ~strcmpi('UTR',M.mut.type));
      else
        fprintf('WARNING: use_nonsilent info not available\n');
        idx = (1:slength(M.mut))';
      end
      if isnumeric(M.mut.gene)
        mutnos = [mutnos, as_row(idx(M.mut.gene(idx)==genes(g)))];
      elseif isfield(M.mut,'gene_idx')
        mutnos = [mutnos, as_row(idx(M.mut.gene_idx(idx)==genes(g)))];
      else
        error('waah!');
      end
     % (added "as_row" 2011-04-06: not sure why it is now necessary --ML
    end
  end
  if ~isempty(muttally), muttally(end-1:end)=[]; end  % remove trailing comma
  gset.muttally{i} = muttally;
  gset.mutnos{i} = mutnos;
end

Mg.gene = gset;
