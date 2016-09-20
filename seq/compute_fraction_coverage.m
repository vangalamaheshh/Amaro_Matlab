function C = compute_fraction_coverage(C)

if C.ncat==1
  C.fcov = sum(C.cov,3) ./ repmat(C.targ.len,1,C.ns);
else
  C.fcov = sum(C.totcov,3) ./ repmat(C.targ.len,1,C.ns);
end
