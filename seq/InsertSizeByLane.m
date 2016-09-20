function InsertSizeByLane(samples,P)
% InsertSizeByLane(samples,P)
%
% Mike Lawrence 2009-06-24

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'InsertSizeByLane_output_dir_suffix','isz');
P = impose_default_value(P,'InsertSizeByLane_output_file_extension','isz');
P = impose_default_value(P,'InsertSizeByLane_java_classname','InsertSizeByLane');
P = impose_default_value(P,'InsertSizeByLane_code','ISZ');
P = impose_default_value(P,'InsertSizeByLane_max_isz','2000');

try

java_classpath = [...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
%   '/seq/software/picard/current/bin/sam-1.05.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

if P.cancer_sample, tn = {'normal','tumor'}; else tn = {'sample'}; end

% make directories

for s=1:length(samples)
  sample = samples{s};
  for i=1:length(tn)
    BAMFile{s,i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '.bam'];
    if ~exist(BAMFile{s,i},'file'), error('Not found: %s',BAMFile{s,i}); end
    OutDir{s,i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '_' P.InsertSizeByLane_output_dir_suffix];
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
    for i=1:length(tn)
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
        file = [OutDir{s,i} '/chr' num2str(c) '.' P.InsertSizeByLane_output_file_extension];
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
        banner = [short_sample P.InsertSizeByLane_code tn{i}(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' ' P.InsertSizeByLane_java_classname ...
           ' ' BAMFile{s,i} ' ' file ' ' num2str(P.InsertSizeByLane_max_isz) ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
  end,end,end % next chr, T/N, sample

  if ~all_done, bwait(jobs); end
end

if did_nothing, fprintf(['All ' P.InsertSizeByLane_code ' files are already up-to-date!\n']); end

catch me, excuse(me); end
