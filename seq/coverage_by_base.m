function coverage_by_base(samples,params)
% coverage_by_base(samples)
%
% samples can be a single sample name or a cell array of sample names
%    e.g.    gbm/0188/wgs    ov/0725/wgs
%
% generates tumor.cbb and normal.cbb directories
%
% Mike Lawrence
% new version 2009-05-15

if ~iscell(samples), samples = {samples}; end

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'blacklist','none');

try

java_classpath = [...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

% make directories

for s=1:length(samples)
  sample = samples{s};
  for i=1:2
    if i==1, tn='normal';else tn='tumor'; end
    BAMFile{s,i} = ['/xchip/cga1/lawrence/' sample '/' tn '.bam'];
    OutDir{s,i} = ['/xchip/cga1/lawrence/' sample '/' tn '_cbb'];
    if ~exist(OutDir{s,i},'dir')
      fprintf('Creating directory %s\n',OutDir{s,i});
      mkdir(OutDir{s,i});
end,end,end

did_nothing = true;
all_done = false;
while(~all_done)
  jobs = [];
  cmds = {}; banners = {};
  all_done = true;
  for s=1:length(samples)
    sample = samples{s};
    short_sample = sample_to_short_sample(sample);
    for i=1:2
      if i==1, tn='normal';else tn='tumor'; end
      if exist(BAMFile{s,i})
        dbam = dir(BAMFile{s,i});
        if isempty(dbam.bytes)
          fprintf('ERROR: %s is broken link!\n',BAMFile{s,i});
          continue
        end
      else
        fprintf('ERROR: %s not found!\n',BAMFile{s,i});
        continue
      end
      for c=1:24
        file = [OutDir{s,i} '/chr' num2str(c) '.cbb'];
        if exist(file,'file')
          dcbb = dir(file);
          if dcbb.bytes>50000000
            if dcbb.datenum>=dbam.datenum
              continue  % file is already up-to-date
            else
              fprintf('%s needs to be refreshed:\n',file);
            end
          else
            fprintf('%s was only partially completed:\n',file);
          end
        end
        did_nothing = false;
        all_done = false;
        banners{end+1,1} = [short_sample 'CBB' tn(1) num2str(c)];
        cmds{end+1,1} = ['"java -Xmx8g -classpath ' java_classpath ' '...
           'CoverageByBase ' BAMFile{s,i} ' ' params.blacklist ' ' ...
           file ' ' num2str(c) '"'];
        % submit jobs as they accumulate
        if length(cmds)>=60
          jobs = [jobs; bsub(cmds,banners)];
          cmds = {}; banners = {};
        end
  end,end,end % next chr, T/N, sample
  if ~isempty(cmds)
    jobs = [jobs; bsub(cmds,banners)];
    cmds = {}; banners = {};
  end
  if ~all_done, bwait(jobs); end
end

if did_nothing, fprintf('All CBB files are already up-to-date!\n'); end

catch me, excuse(me); end
