function Nnc = collapse_genes_to_coding_strand(Nngc,G)
% Nnc = collapse_genes_to_coding_strand(Nngc,G)
%
% Given an array Nngc, with 64 rows, 5 columns (N, n->A, n->C, n->G, n->T), and 1 page per gene
%   and a struct G telling which strand each gene is on,
% Computes array Nnc, with 64 rows and 5 columns.
% Collapses (+)strand genes together with reverse-complemented (-)strand genes.

if size(Nngc,1)~=64, error('size(Nngc,1)~=64'); end
if size(Nngc,2)~=5,  error('size(Nngc,1)~=64'); end
if size(Nngc,3)~=slength(G), error('size(Nngc,3)~=slength(G)'); end

demand_field(G,'strand');

plusstrand = strcmp(G.strand,'+');
minusstrand = strcmp(G.strand,'-');

Nnc = nan(64,5);
compbase('ACGT') = 'TGCA';
X = generate_categ_context65_names();
for i=1:64
  oldname = X.name{i};
  newname = [compbase(oldname(1)) ' in ' compbase(oldname(end)) '_' compbase(oldname(end-2))];
  j = find(strcmp(oldname,X.name));
  k = find(strcmp(newname,X.name));
  Nnc(i,1) = sum(Nngc(j,1,plusstrand),3) + sum(Nngc(k,1,minusstrand),3);
  Nnc(i,2) = sum(Nngc(j,2,plusstrand),3) + sum(Nngc(k,5,minusstrand),3);
  Nnc(i,3) = sum(Nngc(j,3,plusstrand),3) + sum(Nngc(k,4,minusstrand),3);
  Nnc(i,4) = sum(Nngc(j,4,plusstrand),3) + sum(Nngc(k,3,minusstrand),3);
  Nnc(i,5) = sum(Nngc(j,5,plusstrand),3) + sum(Nngc(k,2,minusstrand),3);
end




