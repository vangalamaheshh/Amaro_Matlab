function annotated_genes = annotate_gene_list(genes,gene_lists,markers)
  
  annotations = zeros(length(genes),length(gene_lists));
  
  for j=1:length(gene_lists)
    [M m1 m2] = match_string_sets_hash(genes,gene_lists{j});
    annotations(m1,j) = 1;
  end
  
  annotated_genes = genes;
  for i=1:size(annotations,1)
    for j=1:size(annotations,2)
      if annotations(i,j)
        annotated_genes{i} = cat(2,annotated_genes{i},markers{j});
      end
    end
  end
