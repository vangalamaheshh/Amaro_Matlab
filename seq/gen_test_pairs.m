function X = gen_test_pairs(chr1,start1,end1,strand1,chr2,start2,end2,strand2,count,readlen,insertmean,insertstd,rgrp)
% region1 and region2 should be nearly identical sequence (e.g. from genomicSuperDups.txt)

X = [];
X.rgrp = repmat(rgrp,count,1);
X.num = (1:count)';
z = nan(count,1);
X.chr1 = repmat(chr1,count,1); X.start1 = z; X.end1 = z; X.strand1 = repmat(strand1,count,1);
X.qual1 = round(100*rand(count,1));
X.chr2 = repmat(chr2,count,1); X.start2 = z; X.end2 = z; X.strand2 = repmat(1-strand2,count,1);
X.qual2 = round(100*rand(count,1));
X.flip = round(1*rand(count,1));
X.seq1 = repmat({''},count,1); X.seq2 = X.seq1;

rng = min((end1-start1),(end2-start2)) - readlen - insertmean - 3*insertstd;

for i=1:count
  pos1 = round(rng*rand);
  pos2 = round(pos1 + insertmean + insertstd*randn);
  if strand1==0, X.start1(i) = start1+pos1; else X.start1(i)=end1-pos1-readlen; end
  if strand2==0, X.start2(i) = start2+pos2; else X.start2(i)=end2-pos2-readlen; end
  X.end1(i) = X.start1(i)+readlen;
  X.end2(i) = X.start2(i)+readlen;
  s1 = upper(genome_region(chr1,X.start1(i),X.end1(i)));
  s2 = upper(genome_region(chr2,X.start2(i),X.end2(i)));
  if strand1==0, X.seq1{i} = s1; else X.seq1{i} = rc(s1); end
  if strand2==1, X.seq2{i} = s2; else X.seq2{i} = rc(s2); end
end  
