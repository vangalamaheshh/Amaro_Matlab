function link_bam_maf_wig(samples,fhdirs)
% combines fh2law and link_maf_and_wiggle

if ~iscell(samples), samples = {samples}; end

bd1 = '/xchip/cga1/firehose_output/Individual';
bd2 = '/xchip/cga1/lawrence';

cmd = '';
for i=1:length(samples)
  if ~exist('fhdirs','var')
    dir1 = get_firehose_dir(samples{i});
  else
    if fhdirs{i}(1)=='/'
      dir1 = fhdirs{i};
    else
      dir1 = [bd1 '/' fhdirs{i}];
    end
  end
%  name2 = lower(regexprep(samples{i},'-','/'));
  name2 = samples{i};
  dir2 = [bd2 '/' name2];
  if ~exist(dir2,'dir')
    fprintf('Creating: %s\n',dir2);
    mkdir(dir2);
  end
  % coding wig
  mutdir = dir1;
  if strcmp(dir1(end-2:end),'wgs'), mutdir = [mutdir '/coding']; end
  mutdir = [mutdir '/mut'];
  d = dir([mutdir '/*.wig.txt']);
  if isempty(d)
    fprintf('      No wig file for %s\n', name2);
    keyboard
    w_cmd = [];
  else
    srcfile = [mutdir '/' d(1).name];
    destfile = [dir2 '/somatic_coverage.wig'];
    w_cmd = ['ln -sf ' srcfile ' ' destfile];
  end
  % maf
  d = dir([mutdir '/*.maf.annotated']);
  if isempty(d)
    fprintf('      No maf file for %s\n', name2);
    keyboard
    m_cmd = [];
  else
    srcfile = [mutdir '/' d(1).name];
    destfile = [dir2 '/somatic_mutations.maf'];
    m_cmd = ['ln -sf ' srcfile ' ' destfile];
  end
  % whole-genome 
  if strcmp(dir1(end-2:end),'wgs')
    mutdir = [dir1 '/mut'];
    % whole-genome wig
    d = dir([mutdir '/*.wig.txt']);
    if isempty(d)
      fprintf('      No genomic wig file for %s\n', name2);
      keyboard
      w2_cmd = [];
    else
      srcfile = [mutdir '/' d(1).name];
      destfile = [dir2 '/somatic_coverage_genome.wig'];
      w2_cmd = ['ln -sf ' srcfile ' ' destfile];
    end
    % whole-genome maf
    d = dir([mutdir '/*.maf.annotated']);
    if isempty(d)
      fprintf('      No genomic maf file for %s\n', name2);
      keyboard
      m2_cmd = [];
    else
      srcfile = [mutdir '/' d(1).name];
      destfile = [dir2 '/somatic_mutations_genome.maf'];
      m2_cmd = ['ln -sf ' srcfile ' ' destfile];
    end
  end
  % make links
  fprintf('Linking to bams+maf+wig for %s\n',name2);
  cmd = [cmd ';'...
    'ln -sf ' dir1 '/bam/tumor.bam ' dir2 ';'...
    'ln -sf ' dir1 '/bam/tumor.bai ' dir2 ';'...
    'ln -sf ' dir1 '/bam/normal.bam ' dir2 ';'...
    'ln -sf ' dir1 '/bam/normal.bai ' dir2 ';'...
    'ln -sf ' dir2 '/tumor.bai ' dir2 '/tumor.bam.bai' ';'...
    'ln -sf ' dir2 '/normal.bai ' dir2 '/normal.bam.bai' ';'...
    w_cmd ';' ...
    w2_cmd ';' ...
    m_cmd ';' ...
    m2_cmd ';' ...
  ];
  if length(cmd)>7000, system(cmd); cmd=''; end
end
if ~isempty(cmd), system(cmd); cmd=''; end
