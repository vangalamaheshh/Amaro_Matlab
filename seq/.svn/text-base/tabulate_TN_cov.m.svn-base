function C = tabulate_TN_cov(pat,fsuf,outname)
% calls TabulateTNGrid.java

for i=1:slength(pat)
  dr = ['/xchip/tcga_scratch/lawrence/' pat.dir{i}];
  tf{i} = [dr '/tumor.' fsuf];
  nf{i} = [dr '/normal.' fsuf];
  demand_file(tf{i});
  demand_file(nf{i});
end

java_classpath = ['/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq' ];

cmd={}; banner={};
for i=1:slength(pat)
    outf = ['/xchip/tcga_scratch/lawrence/' pat.dir{i} '/' outname];
    if ~exist(outf,'file')
      cmd{end+1} = ['"java -classpath ' java_classpath ' TabulateTNGrid ' tf{i} ' ' nf{i} ' 127 ' outf '"'];
      banner{end+1} = ['TG' pat.short{i}];
      if length(cmd)>30, bsub(cmd,banner); cmd={}; banner={}; end
end,end,end
if ~isempty(cmd), bsub(cmd,banner); cmd={}; banner={}; end

