function link_maf_and_wiggle(samples,fhdirs)

if ~iscell(samples), samples = {samples}; end

for i=1:length(samples)
  fprintf('Copying MAF and wiggle links for %s\n',samples{i});
  if ~exist('fhdirs','var')
    srcdir = get_firehose_dir(samples{i});
  else
    srcdir = fhdirs{i};
  end
  destdir = ['/xchip/cga1/lawrence/' samples{i}];
  if strcmp(srcdir(end-2:end),'wgs'), srcdir = [srcdir '/coding']; end
  srcdir = [srcdir '/mut'];
  d = dir([srcdir '/*.wig.txt']);
  if isempty(d)
    fprintf('      No wig file for %s\n',samples{i});
    w_cmd = [];
  else
    srcfile = [srcdir '/' d(1).name];
    destfile = [destdir '/somatic_coverage.wig'];
    w_cmd = ['ln -sf ' srcfile ' ' destfile];
  end
  d = dir([srcdir '/*.maf.annotated']);
  if isempty(d)
    fprintf('      No maf file for %s\n',samples{i});
    m_cmd = [];
  else
    srcfile = [srcdir '/' d(1).name];
    destfile = [destdir '/somatic_mutations.maf'];
    m_cmd = ['ln -sf ' srcfile ' ' destfile];
  end
  system([w_cmd ';' m_cmd]);
end

