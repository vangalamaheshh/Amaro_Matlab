function W = loadCoverageFromFWB(fwbname)

import org.broadinstitute.cga.tools.seq.FixedWidthBinary;

fprintf('Loading FWB...\n');
f = FixedWidthBinary(fwbname);
W = cell(24,1);
chrlen = load_chrlen;
for i=1:24, fprintf('chr%d ',i);
  W{i} = f.get(i,1,chrlen(i)*1.1);
  W{i}(W{i}==-1)=0;
end, fprintf('\n');
f.close();

