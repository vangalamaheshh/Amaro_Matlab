function GCContentByLane(samples,P)
% GCContentByLane(samples,P)
%
% Mike Lawrence 2009-07-30

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'GCContentByLane_output_dir_suffix','gc');
P = impose_default_value(P,'GCContentByLane_file_extension','gc');
P = impose_default_value(P,'GCContentByLane_java_classname','GCContentByLane');
P = impose_default_value(P,'GCContentByLane_code','GC');
P = impose_default_value(P,'chrfiles_dir','/xchip/tcga/gbm/analysis/lawrence/genome/hg18');

try

java_classpath = [...
   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

% make directories

for s=1:length(samples)
  sample = samples{s};
  for i=1:2
    if i==1, tn='normal';else tn='tumor'; end
    BAMFile{s,i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '.bam'];
    BAIFile{s,i} = find_bai(BAMFile{s,i});
    OutDir{s,i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_' P.GCContentByLane_output_dir_suffix];
    if ~exist(OutDir{s,i},'dir')
      fprintf('Creating directory %s\n',OutDir{s,i});
      mkdir(OutDir{s,i});
end,end,end

did_nothing = true;
all_done = false;
while(~all_done)
  jobs = [];
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
        file = [OutDir{s,i} '/chr' num2str(c) '.' P.GCContentByLane_file_extension];
        if exist(file,'file')
          dcbb = dir(file);
          if dcbb.bytes>0
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
        banner = [short_sample P.GCContentByLane_code tn(1) num2str(c)];
        if c<23, chrfile = ['chr' num2str(c) '.txt'];
        elseif c==23, chrfile = 'chrX.txt';
        elseif c==24, chrfile = 'chrY.txt';
        else error('invalid chr'); end
        chrfile = fullfile(P.chrfiles_dir,chrfile);
        cmd = ['"java -Xmx2g -classpath ' java_classpath ' ' P.GCContentByLane_java_classname ...
           ' ' BAMFile{s,i} ' ' BAIFile{s,i} ' ' file ' ' chrfile  ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
  end,end,end % next chr, T/N, sample

  if ~all_done, bwait(jobs); end
end

if did_nothing, fprintf(['All ' P.GCContentByLane_code ' files are already up-to-date!\n']); end

catch me, excuse(me); end
