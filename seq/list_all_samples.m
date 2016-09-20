function samples = list_all_samples

basedir = '/xchip/tcga_scratch/lawrence';
samples = {};
d1 = dir(basedir);
for i1=1:length(d1)
  if ~d1(i1).isdir, continue; end
  dir1 = d1(i1).name;   % e.g. gbm
  d2 = dir([basedir '/' dir1]);
  for i2=1:length(d2)
    if ~d2(i2).isdir, continue; end
    dir2 = d2(i2).name;    % e.g. 0188
    d3 = dir([basedir '/' dir1 '/' dir2]);
    for i3=1:length(d3)
      if ~d3(i3).isdir, continue; end
      dir3 = d3(i3).name;    % e.g. wgs
      sample = [dir1 '/' dir2 '/' dir3];
      tbam = dir([basedir '/' sample '/tumor.bam']);
      nbam = dir([basedir '/' sample '/normal.bam']);
      if ~isempty(tbam) && ~isempty(nbam)
        samples = [samples; sample];
end,end,end,end

