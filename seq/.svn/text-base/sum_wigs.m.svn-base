function sum_wigs(samples,outwig)
% Mike Lawrence 2010

javaclasspath('/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq');
import java.io.*;
import java.lang.*;

fprintf('Allocating wiggle room\n');
g = Genome();

for i=1:length(samples)
  wig = ['/xchip/tcga_scratch/lawrence/' samples{i} '/somatic_coverage.wig'];
  fprintf('Loading %d/%d: %s\n',i,length(samples),wig);
  g.loadWiggle(String(wig),String('sum'));
end

fprintf('Saving %s\n',outwig);
g.saveWiggle(outwig);

