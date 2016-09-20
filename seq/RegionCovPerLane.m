function RegionCovPerLane(samples,targfile,outstem,P)
% RegionCovPerLane(samples,targfile,outstem)
%
% samples can be a single sample name or a cell array of sample names
%
% targfile points to a file with one line per target region:  <targname/genename>   <chr>   <start>   <end>
%      must be tab-delimited
%
% outstem specifies the output files that will be created in the sample's home directory:
%    one for tumor, and one for normal
% 
% output file has one line per target,
% breaking down by lane the number of reads overlapping the target:
%     <targname/genename> <chr> <start> <end> <lane0_count> <lane1_count> ... <laneN_count>
%
% Mike Lawrence 2009-06-30
%      modified 2009-10-15

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Last parameter should be a P struct'); end

if ~iscell(samples), samples = {samples}; end

P = impose_default_value(P,'use_LSF',true);

java_classpath = [...
%   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
%   '/seq/software/picard/current/bin/sam-1.05.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

try

did_nothing = true;
all_done = false;
done = false(length(samples),2);
while (~all_done)
 all_done = true;
 jobs = {};
 for s = 1:length(samples)
  sample = samples{s};
  short_sample = sample_to_short_sample(sample);
  direc = ['/xchip/tcga_scratch/lawrence/' sample];
  for i=1:2
    if i==1, tn='tumor'; else tn='normal'; end
    bamfile = [direc '/' tn '.bam'];
    if exist(bamfile,'file')
      dbam = dir(bamfile);
      if isempty(dbam.bytes)
        fprintf('ERROR: %s is broken link!\n',bamfile);
        continue
      end
    else
      fprintf('ERROR: %s not found!\n',bamfile);
      continue
    end
    outfile = [direc '/' outstem '_' tn '.txt'];
    if exist(outfile,'file')
      dout = dir(outfile);
      if dout.datenum >= dbam.datenum
        done(s,i) = true;
        continue;   % final file is already up-to-date
      else
        fprintf('%s needs to be refreshed:\n',outfile);
      end
    end
    % submit per-chromosome jobs
    for c=1:24
      tmpdir = [outfile '_tmp'];
      if ~exist(tmpdir), mkdir(tmpdir); end
      chrfile = [tmpdir '/chr' num2str(c) '.txt'];
      if exist(chrfile,'file')
        dchr = dir(chrfile);
        if dchr.datenum >= dbam.datenum
          continue;   % chrfile is already up-to-date
        end
      end
      did_nothing = false;
      all_done = false;
      banner = [short_sample 'RCL' tn(1) num2str(c)];
      cmd = ['java -classpath ' java_classpath ' '...
         'RegionCovPerLane ' bamfile ' ' targfile ' ' chrfile ' ' num2str(c)];
      if P.use_LSF
        jobs = [jobs;bsub(['"' cmd '"'],banner)];
      else
        system(cmd);
      end
    end % next chromosome
  end % next T/N
 end % next sample
 if (all_done) break; end
 bwait(jobs);
end

%if did_nothing
%  fprintf('All files already up-to-date.\n');
%  return;
%end

% concatenate

for s = 1:length(samples)
  sample = samples{s};
  short_sample = sample_to_short_sample(sample);
  direc = ['/xchip/tcga_scratch/lawrence/' sample];
  for i=1:2
    if done(s,i), continue; end
    if i==1, tn='tumor'; else tn='normal'; end
    fprintf('Concatenating: %s %s\n',sample,tn);
    outfile = [direc '/' outstem '_' tn '.txt'];
    tmpdir = [outfile '_tmp'];
    catfile = [direc '/X' outstem '_' tn '.txt'];
    system(['cat ' tmpdir '/chr*.txt > ' catfile]);
    if exist(catfile,'file')
      system(['mv ' catfile ' ' outfile]);
    else
      fprintf('Problem concatenating files for sample %s %s\n',sample,tn);
      keyboard
    end
  end
end

fprintf('All files complete.\n');

catch me, excuse(me); end
