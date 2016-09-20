function r = openBAMFile(sample,tn)
%  openBAMFile(sample,tn)
%
% sample
%    e.g.    gbm/0188    ov/0725
%
% tn
%    t(umor) or n(ormal)
% -->
%   BAMFile = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '.bam'];



import net.sf.samtools.*;
import java.io.*;
import java.lang.*;

if ~contains(sample,'bam')
  switch lower(tn(1))
   case 't', tn = 'tumor';
   case 'n', tn = 'normal';
   case 's', tn = 'sample';
   otherwise error('tn must be "tumor" or "normal" or "sample"');
  end

  BAMFile = ['/xchip/cga1/lawrence/' sample '/' tn '.bam'];
else
  BAMFile = sample;
end

BAIFile = find_bai(BAMFile);
r = SAMFileReader(File(BAMFile),File(BAIFile));

