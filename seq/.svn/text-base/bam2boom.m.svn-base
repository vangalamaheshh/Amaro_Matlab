function jobs=bam2boom(arg1,arg2,arg3)
% bam2boom(sample[,chrs])
%
% Converts the tumor+normal bams to booms in the conventional directories
%
% or
%
% bam2boom(bamfile,boomdir[,chrs])
%
% Converts the specified bam file to boom format in the specified directory.
% 
% Mike Lawrence 2009-08-27

java_classpath = [...
%   '/seq/software/picard/current/bin/sam-1.05.jar:'...
%  '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
   '/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/analysis_pipeline/tools/classes'];
%   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'];

java_class = 'org.broadinstitute.cga.tools.seq.Bam2Boom';

refdir = '/xchip/tcga/gbm/analysis/lawrence/genome/hg18';

mode = 0;
chrs = [1:24 0];

if ~exist('arg1','var'), error('at least one parameter required'); end
if ~ischar(arg1), error('first parameter must be a string'); end

if exist('arg2','var') && ischar(arg2)
  mode = 1;
  bamfile = arg1;
  boomdir = arg2;
  if exist('arg3','var'), chrs = arg3; end
else
  mode = 2;
  sample = arg1;
  if exist('arg2','var'), chrs = arg2; end
  if exist('arg3','var'), error('check usage'); end
end

if any(chrs<0 | chrs>24), error('chrs should be subset of [1:24 0]'); end

if mode==1   % bam2boom(bamfile,boomdir[,chrs])

  if ~exist(bamfile,'file'), error('%s not found',bamfile); end
  baifile = find_bai(bamfile);
  if ~exist(boomdir,'dir'), fprintf('Creating directory %s\n',boomdir);mkdir(boomdir); end

  ss=[];
  try
    tmp = parse(boomdir,'/([^/]*)/wgs/boom/(tumor|normal)\.boom$',{'smp','tn'});
    if ~isempty(tmp)
      ss = [tmp.smp{1} tmp.tn{1}(1)];
    end
  catch me
  end
  if isempty(ss)
    ss = boomdir;
    ss(1:length(ss)-8)=[];
  end

  dbam = dir(bamfile);

  jobs=[];
  did_nothing = true;
  for ci=1:length(chrs), c=chrs(ci);
    if c>0
      fname = [boomdir '/chr' num2str(c) '.boomindex'];
    else
      fname = [boomdir '/boom.readgroups'];
    end
    if exist(fname,'file')
      dboom = dir(fname);
      if dboom.datenum >= dbam.datenum
        continue;
      end
    end
    banner = [ss 'BB' num2str(c)];
    cmd = ['-E "cd /xchip/cga1" "java -classpath ' java_classpath ' ' java_class ' '...
         bamfile ' ' refdir ' ' boomdir ' ' num2str(c) '"'];
    jobs=[jobs;bsub(cmd,banner)];
    did_nothing = false;
  end

elseif mode==2   % bam2boom(sample[,chrs])

  basedir = '/xchip/cga1/lawrence';

  ss = sample_to_short_sample(sample);

  tn = {'tumor','normal'};

  jobs=[];
  for i=1:length(tn)
    bamfile = [basedir '/' sample '/' tn{i} '.bam'];
    baifile = find_bai(bamfile);
    boomdir = [basedir '/' sample '/' tn{i} '.boom'];
    if ~exist(boomdir,'dir'), mkdir(boomdir); end

    dbam = dir(bamfile);

    did_nothing = true;
    for ci=1:length(chrs), c=chrs(ci);
      if c>0
        fname = [boomdir '/chr' num2str(c) '.boomindex'];
      else
        fname = [boomdir '/boom.readgroups'];
      end
      if exist(fname,'file')
      dboom = dir(fname);
      if dboom.datenum >= dbam.datenum
        continue;
      end
    end
    banner = [ss 'BB' tn{i}(1) num2str(c)];
    cmd = ['"java -classpath ' java_classpath ' ' java_class ' '...
         bamfile ' ' refdir ' ' boomdir ' ' num2str(c) '"'];
    jobs=[jobs;bsub(cmd,banner)];
    did_nothing = false;
   end
  end

else
  error('Unknown mode: check usage');
end

if did_nothing
  fprintf('Boom files already up-to-date\n');
else
  fprintf('All jobs submitted.  40x coverage BAM file takes ~1 hr.\n');
end
