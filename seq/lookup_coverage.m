function [t n] = lookup_coverage(X,fsuf)

if ~isnumeric(X.chr), X.chr=convert_chr(X.chr); end
if ~isnumeric(X.start), X.start = str2double(X.start); end

t = nan(slength(X),1);
n = nan(slength(X),1);

[u ui uj] = unique(X.sample_dir);
for i=1:length(u), fprintf('%d/%d %s\n',i,length(u),u{i});
  d = ['/xchip/tcga_scratch/lawrence/' u{i}];
  idx = find(uj==i);
  c = X.chr(idx);
  p = X.start(idx);
  w = Wiggle([d '/tumor.' fsuf]); t(idx) = w.extract(c,p); w.close();
  w = Wiggle([d '/normal.' fsuf]); n(idx) = w.extract(c,p); w.close();
end, fprintf('\n');

