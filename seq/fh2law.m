function fh2law(samples,fhdirs)

if ~iscell(samples), samples = {samples}; end

bd1 = '/xchip/cga1/firehose_output/Individual';
bd2 = '/xchip/cga1/lawrence';

for i=1:length(samples)
  if ~exist('fhdirs','var')
    dir1 = get_firehose_dir(samples{i});
  else
    dir1 = fhdirs{i};
  end
%  name2 = lower(regexprep(samples{i},'-','/'));
  name2 = samples{i};
  dir2 = [bd2 '/' name2];
  if ~exist(dir2,'dir')
    fprintf('Creating: %s\n',dir2);
    mkdir(dir2);
  end
  fprintf('Linking to BAMs for %s\n',name2);
  system(['ln -sf ' dir1 '/bam/tumor.bam ' dir2 ';'...
  'ln -sf ' dir1 '/bam/tumor.bai ' dir2 ';'...
  'ln -sf ' dir1 '/bam/normal.bam ' dir2 ';'...
  'ln -sf ' dir1 '/bam/normal.bai ' dir2 ';'...
  'ln -sf ' dir2 '/tumor.bai ' dir2 '/tumor.bam.bai' ';'...
  'ln -sf ' dir2 '/normal.bai ' dir2 '/normal.bam.bai']);
end

