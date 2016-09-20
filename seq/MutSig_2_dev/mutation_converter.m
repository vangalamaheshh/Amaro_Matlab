function cDNA_pos = mutation_converter(exons,genomic_pos)

total_matrix = zeros(size(exons,1) + 2, size(exons,1) + 2);
length_column = zeros(size(exons,1), 1);
cumsum_column = zeros(size(exons,1), 1);
cDNA_pos = nan(size(genomic_pos, 1),1); 
for i=1:size(exons,1) 
  length_column(i) = exons(i, 2) - exons(i, 1) + 1;  
end
cumsum_column = cumsum(length_column);
total_matrix = [exons length_column cumsum_column];
for i = 1:size(genomic_pos,1) 
  for j = 1:size(exons,1) 
     if genomic_pos(i)>=total_matrix(j,1) & genomic_pos(i)<=total_matrix(j,2)
      offset = total_matrix(j,2) - genomic_pos(i);
      cDNA_pos(i) = total_matrix(j,4) - offset;
      break; 
    end
  end
end

%%%  if any mutations failed to map to an exon, replace them with NaN
%cDNA_pos(cDNA_pos<1 | cDNA_pos>length) = nan;
