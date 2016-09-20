function W = loadCoverageFromWiggle(wigname)

import org.broadinstitute.cga.tools.seq.Genome;

fprintf('Allocating genome space\n');
g = Genome(true);
fprintf('Loading wiggle file\n');
g.loadWiggle(wigname);
fprintf('\nConverting to Matlab\n');
W = cell(24,1);
chrlen = load_chrlen;
for i=1:24, fprintf('chr%d ',i);
  W{i} = g.getContents(i,1,chrlen(i));
  W{i}(W{i}==-128)=0;
end, fprintf('\n');

