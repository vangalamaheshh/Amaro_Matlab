function dRangerPreprocess(sample,P)
% dRangerPreprocess(sample,P)
%
% (1) Runs all FindWeird jobs
% (2) Runs cat+sort+join on FindWeird jobs to produce "all.weird.pairs"
% (3) Loads "all.weird.pairs" into memory
% (4) Separates and saves two structs as matfiles:
%     "dRanger_input.mat" = all numeric columns
%     "dRanger_seqs.mat" = all string columns (if P.write_weird_pair_sequences is true)
%
% Mike Lawrence 2009-07-23
%      modified 2009-10-16

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Multiple samples not supported'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'blacklist','none');
P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'tumor_only',false);
P = impose_default_value(P,'normal_only',false);
P = impose_default_value(P,'tumor_weird_pair_cutoff_strand_normal',600);
P = impose_default_value(P,'tumor_weird_pair_cutoff_strand_weird',200);
P = impose_default_value(P,'normal_weird_pair_cutoff_strand_normal',450);
P = impose_default_value(P,'normal_weird_pair_cutoff_strand_weird',150);
P = impose_default_value(P,'sample_weird_pair_cutoff_strand_normal',600);
P = impose_default_value(P,'sample_weird_pair_cutoff_strand_weird',200);
P = impose_default_value(P,'write_weird_pair_sequences',false);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR1');
P = impose_default_value(P,'dRangerPreprocess_reads_file','all.weird.reads');
P = impose_default_value(P,'dRangerPreprocess_pairs_file','all.weird.pairs');
P = impose_default_value(P,'dRanger_input_matfile','dRanger_input.mat');
P = impose_default_value(P,'dRanger_seqs_matfile','dRanger_seqs.mat');

java_classpath = [...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sam.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

sort_and_join_script = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sort_and_join.csh';

short_sample = sample_to_short_sample(sample);

if P.cancer_sample
  if P.tumor_only & P.normal_only, error('Cannot specify tumor_only AND normal_only');
  elseif P.tumor_only, tn = {'tumor'};
  elseif P.normal_only, tn = {'normal'};
  else tn = {'tumor','normal'};
  end
else
  tn = {'sample'};
end

if P.write_weird_pair_sequences, writeseqs=1; else writeseqs=0; end

for i=1:length(tn)
  BAMFile{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '.bam'];
  if ~exist(BAMFile{i},'file'), error('Cannot find %s',BAMFile{i}); end
  dbam{i} = dir(BAMFile{i});
  OutDir{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '_' P.dRangerPreprocess_output_dir_suffix];
  if ~exist(OutDir{i},'dir'), mkdir(OutDir{i}); end
end

% FINDWEIRD

did_nothing = true;
all_done = false;
while(~all_done)
  all_done = true;
  jobs = [];
  for c=1:24
    for i=1:length(tn)
      cutoffSN = getfield(P,[tn{i} '_weird_pair_cutoff_strand_normal']);
      cutoffSW = getfield(P,[tn{i} '_weird_pair_cutoff_strand_weird']);
      outfile = [OutDir{i} '/chr' num2str(c) '.weird'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam{i}.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end    
      if ~uptodate
        all_done = false;
        did_nothing = false;
        banner = [short_sample 'FW' tn{i}(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' '...
               'FindWeird ' BAMFile{i} ' ' P.blacklist ' ' num2str(cutoffSW) ' ' num2str(cutoffSN) ...
               ' ' num2str(writeseqs) ' ' outfile ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end

  if ~all_done
    fprintf('Waiting for FindWeird jobs to finish\n');
    bwait(jobs);
  else
    break
  end

end % while(~all_done)

% CONCATENATE

jobs = [];
all_done = false;
while(~all_done)
  all_done = true;
  for i=1:length(tn)
    outfile = [OutDir{i} '/' P.dRangerPreprocess_reads_file];      % all.weird.reads
    dout = dir(outfile);
    uptodate = false;
    if ~isempty(dout)
      if dout.datenum>=dbam{i}.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',outfile); end
    end
    if ~uptodate
      did_nothing = false;
      all_done = false;
      banner = [short_sample 'CAT' tn{i}(1)];
      cmd = ['"cat ' OutDir{i} '/chr*.weird > ' outfile '"'];
      jobs = [jobs;bsub(cmd,banner)];
    end
  end
  if ~all_done
    fprintf('Waiting for concatenate to finish\n');
    bwait(jobs);
  end
end

% SORT+JOIN

jobs = [];
all_done = false;
while(~all_done)
  all_done = true;
  for i=1:length(tn)
    infile = [OutDir{i} '/' P.dRangerPreprocess_reads_file];      % all.weird.reads
    outfile = [OutDir{i} '/' P.dRangerPreprocess_pairs_file];     % all.weird.pairs
    dout = dir(outfile);
    uptodate = false;
    if ~isempty(dout)
      if dout.datenum>=dbam{i}.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',outfile); end
    end
    if ~uptodate
      did_nothing = false;
      all_done = false;
      banner = [short_sample 'STJN' tn{i}(1)];
      cmd = ['-R "rusage[mem=8]" "' sort_and_join_script ' ' infile ' ' outfile '"'];
      jobs = [jobs;bsub(cmd,banner)];
    end
  end
  if ~all_done
    fprintf('Waiting for sort+join to finish\n');
    bwait(jobs);
  end
end

% CONVERT TO STRUCT MATFILES

for i=1:length(tn)
  infile = [OutDir{i} '/' P.dRangerPreprocess_pairs_file];     % all.weird.pairs
  outfile1 = [OutDir{i} '/' P.dRanger_input_matfile];          % dRanger_input.mat
  outfile2 = [OutDir{i} '/' P.dRanger_seqs_matfile];           % dRanger_seqs.mat
  din = dir(infile);
  if isempty(din), error('Not found: %s',infile); end
  dout1 = dir(outfile1);
  dout2 = dir(outfile2);
  if P.write_weird_pair_sequences
    if ~isempty(dout1) && ~isempty(dout2) && dout1.datenum>=din.datenum && dout2.datenum>=din.datenum
      fprintf('Already exist and up-to-date:\n  %s\n  %s\n',outfile1,outfile2);
      continue;
    end
  else    
    if ~isempty(dout1) && dout1.datenum>=din.datenum
      fprintf('Already exist and up-to-date:\n  %s\n',outfile1);
      continue;
    end
  end
  did_nothing = false;
  fprintf('Loading %s\n',infile);
  flds = {'rgrp','namenumber','chr1','start1','end1','strand1','qual1',...
                      'chr2','start2','end2','strand2','qual2','flip'};
  fmt = repmat('%f',1,13);
  if P.write_weird_pair_sequences
    flds = [flds,'seq1','seq2'];
    fmt = [fmt '%s%s'];
  end
  X = load_struct(infile,fmt,0);
  X = rename_fields(X,colx(1:length(flds)),flds);
  if P.write_weird_pair_sequences
    Xseqs = keep_fields(X,{'seq1','seq2'});
    fprintf('Saving %s\n',outfile2);
    save(outfile2,'Xseqs','-v7.3');
    X = rmfield(X,{'seq1','seq2'});
  end
  fprintf('Saving %s\n',outfile1);
  save(outfile1,'X','-v7.3');
end

if did_nothing
  fprintf('All dRangerPreprocess output files already up-to-date\n');
else
  fprintf('dRangerPreprocess finished.\n');
end
  

