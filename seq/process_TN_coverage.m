function C = process_TN_coverage(C)

% collapse categories
C.tumcov = sum(C.tumcov,3); C.ncat = 1;
C.normcov = sum(C.normcov,3); C.ncat = 1;

% collapse to genes
[C.gene.name ui uj] = unique(C.targ.gene);
C.ng = slength(C.gene); C.gtumcov = nan(C.ng,C.ns,C.ncat); C.gnormcov = C.gtumcov;
for i=1:C.ng
  idx = find(uj==i);
  C.gene.chr(i,1) = C.targ.chr(idx(1));
  C.gene.start(i,1) = min(C.targ.start(idx));
  C.gene.end(i,1) = max(C.targ.end(idx));
  C.gene.len(i,1) = sum(C.targ.len(idx));
  C.gene.gc(i,1) = weighted_mean(C.targ.gc(idx),C.targ.len(idx));
  C.gtumcov(i,:) = sum(C.tumcov(idx,:),1);
  C.gnormcov(i,:) = sum(C.normcov(idx,:),1);
end

% convert #bases to #reads (approx)

C.gtumcov = round(C.gtumcov / 76);
C.gnormcov = round(C.gnormcov / 76);



