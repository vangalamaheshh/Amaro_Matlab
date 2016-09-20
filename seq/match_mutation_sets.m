function [nm nc ns nc2 gene_name] = match_mutation_sets(M1,M2,P);
%
% M1 should be TCGA or TSP, etc.
%    must have gene.name, mut.gene, mut.chr, mut.start, mut.end, mut.site, and ng
% M2 should be Cosmic, as loaded by load_cosmic_database()
%    must have gene.name, mut.gene, mut.chr, mut.start, mut.end, mut.site
% P is parameters
%    can have match_margin, default=0
%
% returns (for each gene in M1):
%    nm =  number of unique sites in M1 (TCGA/TSP)
%    nc =  number of unique sites in M2 (Cosmic)
%    ns =  number of sites that overlap between the two data sets
%    nc2 = amount of territory to consider (with each site expanded according to match_margin)
%

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'match_margin',0);

nm = zeros(M1.ng,1);
nc = zeros(M1.ng,1);
ns = zeros(M1.ng,1);
nc2 = zeros(M1.ng,1);

B = P.match_margin; % margin of error for matching coordinates

M2.mut.genename = repmat({''},slength(M2.mut),1);
M2.mut.genename(~isnan(M2.mut.gene)) = M2.gene.name(M2.mut.gene(~isnan(M2.mut.gene)));

if iscell(M1.mut.gene)
  M1.mut.gene = listmap(M1.mut.gene,M1.gene.name);
end

if ~isfield(M1,'use_nonsilent')
  M1.use_nonsilent = 1:slength(M1.mut);
end

fprintf('Matching mutations\n');
fprintf('Gene');

gene_name = M1.gene.name;

for i=1:M1.ng
  if mod(i,1000)==0, fprintf(' %d/%d', i, M1.ng); end
  gname = M1.gene.name{i};
  gno = i;
  m = M1.use_nonsilent(find(M1.mut.gene(M1.use_nonsilent)==gno));
  [tmp x tmp] = unique(M1.mut.site(m));
  m=m(x);
  nm(i) = length(m);
  c = find(strcmp(M2.mut.genename,gname) & ~isnan(M2.mut.chr));
  [tmp x tmp] = unique(M2.mut.site(c));
  c=c(x);
  nc(i) = length(c);

  % compute territory (nc2)
  z = [];
  for k=1:length(c), z = [z,c(k)-B:c(k)+B]; end
  nc2(i) = length(unique(z));

  % find overlaps
  for mi=1:nm(i)
    for ci=1:nc(i)
       if M2.mut.chr(c(ci))==M1.mut.chr(m(mi)) ...
         && M1.mut.start(m(mi)) <= M2.mut.end(c(ci))+B ...
         && M1.mut.end(m(mi)) >= M2.mut.start(c(ci))-B
        ns(i) = ns(i) + 1;
        break;
      end
    end
  end
end
fprintf('\n');
