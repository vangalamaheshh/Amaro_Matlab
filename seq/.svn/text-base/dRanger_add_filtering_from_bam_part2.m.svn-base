function X = dRanger_add_filtering_from_bam_part2(X,tbam,nbam,params)
% Mike Lawrence 2009-09

if ~exist('params','var'), params=[]; end

params=impose_default_value(params,'quiet',1);
params=impose_default_value(params,'refdir','/xchip/tcga/gbm/analysis/lawrence/genome/hg18');
params=impose_default_value(params,'blacklist','none');
params=impose_default_value(params,'skip_end2',false);
params=impose_default_value(params,'maxreads',10000);

import java.io.*;
import java.lang.*;

F = nan(slength(X),16);
col = [1 9;5 13];

bam = {tbam nbam};
chr = [X.chr1 X.chr2];
mn = [X.min1 X.min2];
mx = [X.max1 X.max2];

for b=1:2
  demand_file(bam{b});
  BG{b} = BamGrasp();
  if params.quiet, BG{b}.setQuietModeOn; else BG{b}.setQuietModeOff; end
  BG{b}.openFile(String(bam{b}),String(params.blacklist),String(params.refdir));
  if ~isempty(params.maxreads), BG{b}.set_maxReads(params.maxreads); end
end

rstep = 1;
for i=1:slength(X), if ~mod(i,rstep), fprintf('%d/%d ',i,slength(X)); end
  for b=1:2
    for e=1:2
      R = BamGrasp_load_region(BG{b}, chr(i,e), mn(i,e), mx(i,e), params);
      cov = size(R,1);
      avgnmm = mean(R(:,7));
      fmapqz = mean(R(:,8)==0);
      nuwp = subfunction_nuwp(R,chr(i,e),mn(i,e));
      F(i,col(b,e)+[0:3]) = [cov fmapqz avgnmm nuwp];
    end
  end
end

for b=1:2
  BG{b}.closeFile();
end

X.coverageT1 = F(:,1);
X.coverageN1 = F(:,5);
X.fmapqzT1 = F(:,2);
X.fmapqzN1 = F(:,6);
X.avgnmmT1 = F(:,3);
X.avgnmmN1 = F(:,7);
X.nuwpT1 = F(:,4);
X.nuwpN1 = F(:,8);

if ~params.skip_end2
  X.coverageT2 = F(:,9);
  X.coverageN2 = F(:,13);
  X.fmapqzT2 = F(:,10);
  X.fmapqzN2 = F(:,14);
  X.avgnmmT2 = F(:,11);
  X.avgnmmN2 = F(:,15);
  X.nuwpT2 = F(:,12);
  X.nuwpN2 = F(:,16);
end

return

  function nuwp = subfunction_nuwp(R,chr,pos)
    WEIRDPAIR_THRESHOLD = 10000;
    WEIRDPAIR_DISTINCTNESS_THRESHOLD = 2000;
    wp = R(:,10:12);                                      % chr start strand
    wp = wp(wp(:,1)>0,:);                                 % remove unpaired and unmapped-pair
    nw = find(wp(:,1)==chr & abs(wp(:,2)-pos)<WEIRDPAIR_THRESHOLD);     % remove non-weird
    wp = wp(setdiff(1:size(wp,1),nw),:);
    wp(:,2) = round(wp(:,2) / WEIRDPAIR_DISTINCTNESS_THRESHOLD);    % smooth to windows of 2kB
    nuwp = size(unique(wp,'rows'),1);
  end

end   % main function
