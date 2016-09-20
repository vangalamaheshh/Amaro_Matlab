function b = binsum(a,ybin,xbin)
% given matrix a and sets of bins ybin and xbin, reformats into b, with bins of summed counts
%
% Mike Lawrence 2010-02-03

xbin = as_column(xbin);
ybin = as_column(ybin);

xto = nan(size(a,2),1);
for i=1:size(a,2)
  idx = find(i>=xbin,1,'last');
  if isempty(idx), xto(i)=length(xbin); else xto(i)=idx; end
end

yto = nan(size(a,1),1);
for i=1:size(a,1)
  idx = find(i>=ybin,1,'last');
  if isempty(idx), yto(i)=length(xbin); else yto(i)=idx; end
end

s1 = zeros(length(ybin),size(a,2));
for i=1:length(ybin)
  s1(i,:) = sum(a(yto==i,:),1);
end

b = zeros(length(ybin),length(xbin));
for i=1:length(xbin)
  b(:,i) = sum(s1(:,xto==i),2);
end
