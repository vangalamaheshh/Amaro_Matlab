function R = annotate_regions_that_are_Ig_loci(R)
% annotates
% three immunoglobulin loci, with somatic hypermutation in B-cell lineages
%
% from Mike Chapman

require_fields(R,{'chr','start','end'});

ig=[]; ig.chr = [2;14;22]; ig.start = [88500000;105000000;21200000]; ig.end = [95800000;inf;21700000];
ig.name = {'IgK','IgH','IgL'};

R.ig = false(slength(R),1);
for i=1:slength(ig)
  idx = (R.chr==ig.chr(i) & R.start<=ig.end(i) & R.end>=ig.start(i));
  R.ig(idx) = true;
end
