function c = count_mf_coverage(X,chr,st,en)
%
% given multifind results in struct X
%   with fields: chr1,start1,end1,chr2,start2,end2
%
% and target region defined by chr,st,en
%
% returns coverage counts over the target region
%   i.e. number of reads overlapping each position (both pairmates counted)
%
% Mike Lawrence 2009-03-14

c = zeros(en-st+1,1);
for i=1:slength(X)
  if X.chr1(i)==chr
    first = max(1,X.start1(i)-st+1);
    last = min(en-st+1,X.end1(i)-st+1);
    if first>last, first=2; last=1; end
    c(first:last) = c(first:last) + 1;
  end
  if X.chr2(i)==chr
    first = max(1,X.start2(i)-st+1);
    last = min(en-st+1,X.end2(i)-st+1);
    if first>last, first=2; last=1; end
    c(first:last) = c(first:last) + 1;
  end
end

