function C = collapse_capture_coverage_to_genes(C)
[C.gene.name ui uj] = unique(C.targ.gene);
C.ng = slength(C.gene); C.gcov = nan(C.ng,C.ns,C.ncat);
for i=1:C.ng, if ~mod(i,1000), fprintf('%d/%d ',i,C.ng); end
  idx = find(uj==i);
  C.gene.chr(i,1) = C.targ.chr(idx(1));
  C.gene.start(i,1) = min(C.targ.start(idx));
  C.gene.end(i,1) = max(C.targ.end(idx));
  C.gene.len(i,1) = sum(C.targ.len(idx));
  C.gene.gc(i,1) = weighted_mean(C.targ.gc(idx),C.targ.len(idx));
  for c=1:C.ncat, C.gcov(i,:,c) = sum(C.cov(idx,:,c),1); end
end,fprintf('\n');
C.fgcov = sum(C.gcov,3) ./ repmat(C.gene.len,1,C.ns);
