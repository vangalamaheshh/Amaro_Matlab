function print_mutation_stats(M)

fprintf('\nMutation statistics:\n');
pgok = find(~isnan(M.mut.patient) & ~isnan(M.mut.gene));
fprintf('  %d genes, %d samples, %d mutations\n', M.ng, M.np, length(pgok));

fprintf('\nFiltering results:\n');
count(M.mut.filtered(pgok));

fprintf('\nBreakdown by center\n');
% remove weird Baylor notations
center = regexprep(M.mut.center, 'TCG_\d*_\d*._SD3_\d\d','BCM');
count(center(M.use));

fprintf('\nBreakdown by type and evidence\n\n');
type = regexprep(M.mut.type, 'Frameshift_...','Frameshift_Indel');
type = regexprep(type, 'Inframe_...','Inframe_Indel');
val = M.use(grep('Validated', M.mut.evidence(M.use), 1));
ver = M.use(grep('Verified', M.mut.evidence(M.use), 1));
seq = M.use(grep('Sequenced', M.mut.evidence(M.use), 1));
categ = { 'Missense'; 'Synonymous'; 'Nonsense'; 'Frameshift_Indel'; 'Inframe_Indel'; 'Splice_site' };
fprintf('%-16s  %8s  %8s  %8s  %8s\n\n', '','Validated','Verified','Sequenced','Total');
for c=1:length(categ)
  valc = length(find(strcmp(categ{c},type(val))));
  verc = length(find(strcmp(categ{c},type(ver))));
  seqc = length(find(strcmp(categ{c},type(seq))));
  totc = valc+verc+seqc;
  fprintf('%-16s  %8d  %8d  %8d  %8d\n', categ{c}, valc, verc, seqc, totc);
end
fprintf('\n%-16s  %8d  %8d  %8d  %8d\n\n','Total', length(val), length(ver), length(seq), length(M.use));


end
