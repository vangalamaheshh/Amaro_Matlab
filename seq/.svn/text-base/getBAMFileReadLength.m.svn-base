function len = getBAMFileReadLength(sample,P)
%  getBAMFileReadLength(sample)
%

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'cancer_sample',true);

%javaclasspath('/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar')
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;

if P.cancer_sample
  rt = openBAMFile(sample,'t');
  it = rt.iterator;
  for i=1:100,lt(i,1)=it.next.getReadLength;end

  rn = openBAMFile(sample,'n');
  in = rn.iterator;
  for i=1:100,ln(i,1)=it.next.getReadLength;end

  len = nanmedian([lt;ln]);
else
  rs = openBAMFile(sample,'s');
  is = rs.iterator;
  for i=1:100,ls(i,1)=is.next.getReadLength;end

  len = nanmedian(ls);
end

