function E = convert_to_exon_list(G)
% E = convert_to_exon_list(G)
%
% converts genelist G (from get_refseq_introns_and_exons)
% to E = list of exons (coding sequence only)

E = cell(slength(G),1);
for i=1:slength(G)
  E{i}.gene = repmat({G.name{i}},G.n_exons(i),1);
  E{i}.chr = G.chr(i) * ones(G.n_exons(i),1);
  E{i}.start = G.exon_starts{i};
  E{i}.end = G.exon_ends{i};
end
E = concat_structs(E);
