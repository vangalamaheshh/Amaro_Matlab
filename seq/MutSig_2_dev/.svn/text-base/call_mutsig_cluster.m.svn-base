function R = call_mutsig_cluster(M, P)

if ~iscell(M), M={M}; end

if ~exist('P','var'), P=[]; end

if ~isfield(P,'genes_to_process')
  P.genes_to_process = -1;  % process all genes
else
  if ischar(P.genes_to_process)
    P.genes_to_process = { P.genes_to_process };
  end
  if ~iscell(P.genes_to_process)
    error('P.genes_to_process should be a list of gene names');
  end
end

%% failsafe check: make sure target list is exons, not whole transcripts
if mean(M{1}.cov.targ.len)>5000
  error('Looks like target list is whole transcripts: this will cause problems for MutSig2.');
end

% PREPROCESS EACH DATASET

nsets = length(M);
setnames = cell(nsets,1);
for si=1:nsets
    setnames{si} = M{si}.name;
    fprintf('Restricting to nonsilent coding mutations\n');
    if isfield(M{si},'use_nonsilent')
      % keep only use_nonsilent mutations
      M{si}.mut = reorder_struct(M{si}.mut,M{si}.use_nonsilent);
      M{si} = rmfield(M{si},'use_nonsilent');
    elseif isfield(M{si}.mut,'is_silent') && isfield(M{si}.mut,'is_coding')
      M{si}.mut = make_numeric(M{si}.mut,{'is_coding','is_silent'});
      M{si}.mut = reorder_struct(M{si}.mut,M{si}.mut.is_coding & ~M{si}.mut.is_silent);
    elseif isfield(M{si}.mut,'type')
      idx = grepi('missense|nonsense|splice|non.?stop|in.?frame|frame.?shift|de.?novo',M{si}.mut.type);
      M{si}.mut = reorder_struct(M{si}.mut,idx);
    else
      error('Don''t know how to select nonsilent coding mutations!');
    end
end

builds = cell(nsets, 1);
for si=1:nsets, builds{si} = M{si}.build; end

%% make a gene list sorted by total number of mutations 
fprintf('Sorting genes by total number of mutations....\n');
ghist = zeros(M{1}.ng,1);
M{si}.mut.gene = M{si}.mut.gene_idx;
for si=1:nsets, ghist = ghist + histc(M{si}.mut.gene,1:M{1}.ng); end
[tmpzz, gene_order] = sort(ghist,'descend');

temp_structs = cell(nsets,1);
for si = 1:nsets, temp_structs{si} = M{si}.mut; end 
all_muts = concat_structs(temp_structs);

%% Filter sorted gene list for only genes mutated at least twice
fprintf('Processing gene list...\n');
[mutated_genes mgi mgj] = unique(nansub(M{1}.gene.name,all_muts.gene));
mgih = histc(mgj,1:length(mutated_genes));

% output [REPORT] SKIP line for each gene with <2 mutations
skipped_genes = mutated_genes(mgih<2);
for i=1:length(skipped_genes)
  fprintf('[REPORT]  %-11s   SKIP: fewer than two mutations\n',skipped_genes{i});
end

mutated_genes = mutated_genes(mgih>=2);

if iscell(P.genes_to_process)
  mutated_genes = intersect(mutated_genes,P.genes_to_process);
end

idx = listmap(mutated_genes, M{1}.gene.name);
idx = idx(~isnan(idx));
for i = 1:length(gene_order)
  if ~ismember(gene_order(i), idx) 
    gene_order(i) = nan; 
  end
end
gene_order = gene_order(~isnan(gene_order));

% Set all output to NaN

R = [];
R.gene_name = M{1}.gene.name;
z = nan(length(M{1}.gene.name), 1);
R.nperm = z;
R.p_clust = z;
R.p_cons = z;
R.p_joint = z;

% Determine if there's nothing to be done

if isempty(gene_order)
  fprintf('In this run there are no genes with >=2 mutations: nothing to do!\n');

else
  
  % Run the clustering algorithm for each gene and fill R struct with p value fields.

  count = 0; 
  for i = as_row(gene_order)

    gname = M{1}.gene.name{i};
    P.gene = i; 
    P.setnames = setnames;

    try
      Ri = mini_driver(M, P);
      R.nperm(i) = Ri.nperm;
      R.p_clust(i) = Ri.p_clust;
      R.p_cons(i) = Ri.p_cons;
      R.p_joint(i) = Ri.p_joint;
    catch me
      disp(me); disp(me.message);
      fprintf('ERROR in mini_driver with gene %s\n',gname);
      keyboard
    end
    
    count = count + 1;
    percent_done = (count/length(gene_order))*100;
    fprintf('Current progress: %.0f percent done.\n', percent_done);
    
  end 

end



fprintf('Finished call_mutsig_cluster\n');





