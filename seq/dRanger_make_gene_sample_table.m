function X = dRanger_make_gene_sample_table(x,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','hg18');

R = load_refseq(P.build)
G = []; [G.name ui] = unique(R.gene);
G.transcript = R.transcript(ui);
G.chr = R.chr(ui);
G.strand = R.strand(ui);
G.tx_start = R.tx_start(ui);
G.tx_end = R.tx_end(ui);
G.tx_len = R.tx_len(ui);
G.code_len = R.code_len(ui);

y = x;
exclude_igr = true;
if exclude_igr
  idx=grep('^IGR',y.site1,1); y.gene1(idx) = repmat({'---'},length(idx),1);
  idx=grep('^IGR',y.site2,1); y.gene2(idx) = repmat({'---'},length(idx),1);
end
y.gidx1 = listmap(y.gene1,G.name);
y.gidx2 = listmap(y.gene2,G.name);
pat=[];
[pat.name upi upj] = unique(y.individual);
q = false(slength(G),slength(pat));
fprintf('Patient ');
for i=1:slength(pat), fprintf('%d/%d ',i,slength(pat));
  yi = reorder_struct(y,upj==i);
  for j=1:slength(G)
    q(j,i) = any(yi.gidx1==j | yi.gidx2==j);
end,end,fprintf('\n');
sq = sum(q>0,2);
[tmp ord] = sort(sq,'descend');

X = [];
X.pat = pat;
X.gene = reorder_struct(G,ord);
X.has_rearr = q(ord,:);
X.gene.n_samps_rearr = sq(ord);

