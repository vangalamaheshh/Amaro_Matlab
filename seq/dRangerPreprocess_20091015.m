function dRangerPreprocess(sample,P)
% dRangerPreprocess(sample,P)
%
% Mike Lawrence 2009-07-23

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Multiple samples not supported'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'tumor_only',false);
P = impose_default_value(P,'normal_only',false);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
P = impose_default_value(P,'dRangerPreprocess_reads_file','all.weird.reads');
P = impose_default_value(P,'dRangerPreprocess_pairs_file','all.weird.pairs');
P = impose_default_value(P,'dRangerPreprocess_nums_file','all.weird.pairs.nums');
P = impose_default_value(P,'dRangerPreprocess_seqs_file','all.weird.pairs.seqs');


P = impose_default_value(P,'dRangerRun_input_matfile','dRangerRun_input.mat');


P = impose_default_value(P,'coverage_window_size',1000);
P = impose_default_value(P,'coverage_filter_flag',0);


P = impose_default_value(P,'CatSortJoin',true);         % whether to run Cat+SortJoin
P = impose_default_value(P,'SegregateSeqs',true);       % whether to run SegregateSeqs
P = impose_default_value(P,'perform_stringent_realign',false);

try

java_classpath = [...
   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

sort_and_join_script = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/sort_and_join.csh';
segregate_seqs_script = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/segregate_seqs.pl';

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

for i=1:length(tn)
  BAMFile{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '.bam'];
  if ~exist(BAMFile{i},'file'), error('Cannot find %s',BAMFile{i}); end
  dbam{i} = dir(BAMFile{i});
  BAIFile{i} = find_bai(BAMFile{i});
  OutDir{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn{i} '_' P.dRangerPreprocess_output_dir_suffix];
  if ~exist(OutDir{i},'dir'), mkdir(OutDir{i}); end
end

did_nothing = true;
all_done = false;
while(~all_done)
 all_done = true;
 jobs = [];

 % FindWeird
 if P.FindWeird
  for c=1:24
    for i=1:length(tn)
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
        banner = [short_sample 'FW2' tn{i}(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' '...
               'FindWeird2 ' BAMFile{i} ' ' BAIFile{i} ' ' outfile ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end
 end % if P.FindWeird

 % MeasureCoverage (OBSOLETE)
 if P.MeasureCoverage
  for c=1:24
    for i=1:length(tn)
      outfile = [OutDir{i} '/chr' num2str(c) '.cov'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam{i}.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'MCV' tn{i}(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' '...
               'MeasureCoverage ' BAMFile{i} ' ' BAIFile{i} ' ' outfile ' ' ...
                num2str(P.coverage_filter_flag) ' ' num2str(P.coverage_window_size) ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end
 end % if P.MeasureCoverage

 % CountMQZ (OBSOLETE)
 if P.CountMQZ
  for c=1:24
    for i=1:length(tn)
      outfile = [OutDir{i} '/chr' num2str(c) '.mqz'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam{i}.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'MQZ' tn{i}(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' '...
                'CountMQZ ' BAMFile{i} ' ' BAIFile{i} ' ' outfile ' ' ...
                num2str(P.coverage_window_size) ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end
 end % if P.CountMQZ

 % MakeLanelist
 if P.MakeLanelist
  for i=1:length(tn)
    outfile = [OutDir{i} '/lanelist.txt'];
    dout = dir(outfile);
    uptodate = false;
    if ~isempty(dout)
      if dout.datenum>=dbam{i}.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',outfile); end
    end
    if ~uptodate
      all_done = false;
      did_nothing = false;
      banner = [short_sample 'LANES' tn{i}(1)];
      cmd = ['"java -classpath ' java_classpath ' '...
             'MakeLanelist ' BAMFile{i} ' ' outfile '"'];
      jobs = [jobs;bsub(cmd,banner)];
    end
  end
 end % if P.MakeLanelist

 if ~all_done
   fprintf('Waiting for dRangerPreprocess jobs to finish\n');
   bwait(jobs);
 else
   break
 end

end % while(~all_done)


% CatSortJoin
if P.CatSortJoin

  % CONCATENATE
  jobs = [];
  all_done = false;
  while(~all_done)
    all_done = true;
    for i=1:length(tn)
      outfile = [OutDir{i} '/' P.dRangerPreprocess_reads_file];
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
        cmd = ['"' sort_and_join_script ' ' infile ' ' outfile '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
    if ~all_done
      fprintf('Waiting for sort+join to finish\n');
      bwait(jobs);
    end
  end
end % if P.CatSortJoin

if P.SegregateSeqs

  % SEGREGATE SEQS
  jobs = [];
  all_done = false;
  while(~all_done)
    all_done = true;
    for i=1:length(tn)
      infile = [OutDir{i} '/' P.dRangerPreprocess_pairs_file];     % all.weird.pairs
      outfile1 = [OutDir{i} '/' P.dRangerPreprocess_nums_file];    % all.weird.pairs.nums
      outfile2 = [OutDir{i} '/' P.dRangerPreprocess_seqs_file];    % all.weird.pairs.seqs
      dout1 = dir(outfile1); dout2 = dir(outfile2);
      uptodate = false;
      if ~isempty(dout1) & ~isempty(dout2)
        if dout1.datenum>=dbam{i}.datenum && dout2.datenum>=dbam{i}.datenum, uptodate = true;
        else fprintf('%s and/or %s need to be refreshed:\n',outfile1, outfile2); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'SEGR' tn{i}(1)];
        cmd = ['"perl ' segregate_seqs_script ' ' infile ' ' outfile1 ' ' outfile2 '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
    if ~all_done
      fprintf('Waiting for segregate_seqs to finish\n');
      bwait(jobs);
    end
  end
end % if P.SegregateSeqs

if did_nothing
  fprintf('All dRangerPreprocess output files already up-to-date\n');
else
  fprintf('dRangerPreprocess finished.\n');

dRanger_stringent_realign(sample,P);
  
catch me, excuse(me); end
