function R = extract_region_coverage_from_wig(R,covwig)

demand_fields(R,{'chr','start','end'});
nr=slength(R);

g = org.broadinstitute.cga.tools.seq.Genome();
g.loadWiggle(covwig);

R.cov = nan(nr,1);
for i=1:nr, if ~mod(i,100000), fprintf('%d/%d ',i,nr); end
  try
    R.cov(i) = g.getSumContents(R.chr(i),R.start(i),R.end(i));
  catch me
    fprintf('Problem with %d:%d-%d\n',R.chr(i),R.start(i),R.end(i));
  end
end,fprintf('\n');
