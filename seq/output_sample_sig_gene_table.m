function output_sample_sig_gene_table(M,P)
% outputs table with:
%    one row per significantly mutated gene
%    one column per sample
%    cell = number of mutations in gene/sample

if exist('P','var') && ischar(P)
  tmp=P;
  P=[];
  P.sample_siggene_table_filename = tmp;
end

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'sample_siggene_table_filename','*required*');
P = impose_default_value(P,'sample_siggene_table_q_cutoff',0.2);

if ~isfield(M,'n_nonsilent')
  fprintf('no M.n_nonsilent: skipping output_sample_sig_gene_table\n');
  return
end

gidx = find(M.gene.qval<=P.sample_siggene_table_q_cutoff);
[tmp ord] = sort(M.gene.qval(gidx));
gidx = gidx(ord);

out = fopen(P.sample_siggene_table_filename,'wt');
fprintf(out,'gene');
for i=1:M.np, fprintf(out,'\t%s',M.patient.name{i}); end
fprintf(out,'\n');
for i=1:length(gidx)
  fprintf(out,'%s',M.gene.name{gidx(i)});
  for j=1:M.np, fprintf(out,'\t%d',M.n_nonsilent(gidx(i),M.TOT,j)); end
  fprintf(out,'\n');
end
fclose(out);
