function audit_getInferredInsertSize(sample,tn)

f = openBAMFile(sample,tn);
i = f.queryContained('chr1',50000000,60000000);
n = 1;
while 1
  if ~mod(n,100), keyboard; end
  x = i.next;
  if x.getMappingQuality==0, continue; end
  if x.getNotPrimaryAlignmentFlag, continue; end
  if ~x.getReadPairedFlag, continue; end
  if ~x.getProperPairFlag, continue; end
  if x.getMateUnmappedFlag, continue; end
  chr1 = x.getReferenceIndex.doubleValue; st1 = x.getAlignmentStart;
  en1 = x.getAlignmentEnd; str1 = x.getReadNegativeStrandFlag;
  chr2 = x.getMateReferenceIndex.doubleValue; st2 = x.getMateAlignmentStart;
  str2 = x.getMateNegativeStrandFlag;
  isize = x.getInferredInsertSize;
  n=n+1;
  fprintf('%d (%d) %d-%d  //  %d (%d) %d-   ////  %d\n',...
    chr1,str1,st1,en1,chr2,str2,st2,isize);
end
