function T = generate_target_list(P)
% generate_target_list(P)
%
% generates list of targets from refseq
%  -- can include noncoding regions
%  -- can restrict to only conserved regions

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'refseqbuild',P.build);
P = impose_default_value(P,'genomeregionbuild',P.build);
P = impose_default_value(P,'splicesiteflank',2);
P = impose_default_value(P,'include_exons',false);
P = impose_default_value(P,'include_UTRs',true);
P = impose_default_value(P,'include_introns',true);
P = impose_default_value(P,'include_promoters',true);
P = impose_default_value(P,'imputed_promoter_length',3000);
P = impose_default_value(P,'restrict_to_conserved_regions',false);
P = impose_default_value(P,'conservation_matfile',[]);
if P.restrict_to_conserved_regions && isempty(P.conservation_matfile)
  error('must specify P.conservation_matfile');
end

PP = [];
PP.add_promoters = P.include_promoters;
PP.imputed_promoter_length = P.imputed_promoter_length;
G = get_refseq_introns_and_exons(P.refseqbuild,P.splicesiteflank,PP);

if P.restrict_to_conserved_regions
  fprintf('Loading conserved regions from %s\n', P.conservation_matfile);
  C = load(P.conservation_matfile,'R'); C=C.R;
  cidx = cell(24,1);
  for i=1:24; cidx{i} = find(C.chr==i); end
end

T = cell(slength(G),1);
fprintf('Gene: ');
for g=1:slength(G), if ~mod(g,1000), fprintf('%d/%d ',g,slength(G)); end
  % find included regions of this transcript
  chr = G.chr(g); if chr<1 || chr>24, continue; end
  if P.include_UTRs & P.include_introns & P.include_exons
   se = [G.gene_start(g) G.gene_end(g)];
  elseif ~P.include_promoters & ~P.include_UTRs & P.include_introns & P.include_exons
   se = [G.code_start(g) G.code_end(g)];
  elseif ~P.include_promoters & ~P.include_UTRs & ~P.include_introns & P.include_exons
   se = [G.exon_starts{g} G.exon_ends{g}];
  elseif ~P.include_promoters & ~P.include_UTRs & P.include_introns & ~P.include_exons
   se = [G.intron_starts{g} G.intron_ends{g}];
  elseif P.include_UTRs & P.include_introns & ~P.include_exons
   se = [];
   if G.gene_start(g)<G.code_start(g)
     se = [se; G.gene_start(g) G.code_start(g)-1];
   end
     se = [se; G.intron_starts{g} G.intron_ends{g}];
   if G.code_end(g)<G.gene_end(g)
     se = [se; G.code_end(g)+1 G.gene_end(g)];
   end
  else
    error('Unsupported inclusion mode');
  end
  % if specified, restrict to conserved regions
  if P.restrict_to_conserved_regions
    idx = cidx{chr};
    idx = idx(C.end(idx)>=se(1,1) & C.start(idx)<=se(end,2));
    cse = [C.start(idx) C.end(idx)];
    se = interval_list_intersect(se,cse);
  end
  n = size(se,1); if n==0, continue; end
  T{g}.gene = repmat({G.name{g}},n,1);
  T{g}.chr = chr*ones(n,1);
  T{g}.start = se(:,1);
  T{g}.end = se(:,2);
end
fprintf('\n');

T = concat_structs(T);

% add GC content
fprintf('Annotating GC content\n');
T.gc = get_gc_content(T.chr,T.start,T.end,P.build);

T.len = T.end-T.start+1;

T.mem = ones(slength(T),1);
