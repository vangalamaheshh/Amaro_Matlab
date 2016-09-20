function G = load_gene_lengths

build = 'hg18'

R = load_refseq(build)
R.chr = convert_chr(R.chr);
R = reorder_struct(R,~isnan(R.chr));

G = [];
[G.name ui uj] = unique(R.gene);
ng = length(G.name);
for i=1:ng, if ~mod(i,1000), fprintf('%d/%d ',i,ng); end
  idx = find(uj==i);
  G.chr(i,1) = mode(R.chr(idx));
  idx = idx(R.chr(idx)==G.chr(i,1));
  G.str{i,1} = concat(unique(R.strand(idx)),'/');
  G.start(i,1) = min(R.tx_start(idx));
  G.end(i,1) = max(R.tx_end(idx));
  G.code_len(i,1) = max(R.code_len(idx));
end, fprintf('\n');

G.len = G.end-G.start+1;

