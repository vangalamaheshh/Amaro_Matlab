function CoverageByBaseAndStrand(bamfile,baifile,outdir)

while(1)
  jobs = [];

  all_done = true;
  for c=1:24
    file = [outdir '/chr' num2str(c) '.cov'];
    if exist(file,'file')
      d = dir(file);
      if d.bytes>50000000,continue; end % already done
    end
    all_done = false;
    cmd = ['"java -Xmx2g -classpath /xchip/tcga/gbm/analysis/'...
         'lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:' ...
         '/xchip/tcga/gbm/analysis/lawrence/sam CoverageByBaseAndStrand ' ...
         bamfile ' ' baifile ' '...
         file ' ' num2str(c) '"'];
    jobs = [jobs;bsub(cmd)];
  end

  if all_done
    break;
  else
    fprintf('Waiting for CoverageByBase to finish\n');
    bwait(jobs);
  end
end
