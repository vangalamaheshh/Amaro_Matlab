function dRangerPreprocess(sample,P)
% dRangerPreprocess(sample,P)
%
% old parameter style: dRangerPreprocess(sample,coverage_filter_flag,coverage_window_size)
%
% Mike Lawrence 2009

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Multiple samples not supported'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'cancer_sample',true);
P=impose_default_value(P,'coverage_window_size',1000);
P=impose_default_value(P,'coverage_filter_flag',0);
P=impose_default_value(P,'FindWeird',true);           % whether to run FindWeird
P=impose_default_value(P,'MeasureCoverage',true);     % whether to run MeasureCoverage
P=impose_default_value(P,'MakeLanelist',true);        % whether to run MakeLanelist
P=impose_default_value(P,'CountMQZ',true);            % whether to run CountMQZ
P=impose_default_value(P,'CatSortJoin',true);         % whether to run Cat+Sort+Join

%if ~exist('coverage_window_size','var') || isempty(coverage_window_size), coverage_window_size = 1000; end
%if ~exist('coverage_filter_flag','var') || isempty(coverage_filter_flag), coverage_filter_flag = 0; end

try

java_classpath = [...
   '/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar:'...
   '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq'...
];

joinweird_perlscript = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/joinweird.pl';


short_sample = sample_to_short_sample(sample);

for i=1:2
  if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
  else if i==1, tn='sample';else break; end, end
  BAMFile{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '.bam'];
  BAIFile{i} = find_bai(BAMFile{i});
  OutDir{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_dR'];
  if ~exist(OutDir{i},'dir'), mkdir(OutDir{i}); end
end

did_nothing = true;
all_done = false;
while(~all_done)
 all_done = true;
 jobs = [];
 for i=1:2
  if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
  else if i==1, tn='sample';else break; end, end

  dbam = dir(BAMFile{i});

  if P.FindWeird
    % FindWeird
    for c=1:24
      outfile = [OutDir{i} '/chr' num2str(c) '.weird'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        all_done = false;
        did_nothing = false;
        banner = [short_sample 'FWE' tn(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' '...
               'FindWeird ' BAMFile{i} ' ' BAIFile{i} ' ' outfile ' ' num2str(c) ' mqual"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end % if P.FindWeird

  if P.MakeLanelist
    % MakeLanelist
    outfile = [OutDir{i} '/lanelist.txt'];
    dout = dir(outfile);
    uptodate = false;
    if ~isempty(dout)
      if dout.datenum>=dbam.datenum, uptodate = true;
      else fprintf('%s needs to be refreshed:\n',outfile); end
    end
    if ~uptodate
      all_done = false;
      did_nothing = false;
      banner = [short_sample 'LANES' tn(1)];
      cmd = ['"java -classpath ' java_classpath ' '...
             'MakeLanelist ' BAMFile{i} ' ' outfile '"'];
      jobs = [jobs;bsub(cmd,banner)];
    end
  end % if P.MakeLanelist

  if P.MeasureCoverage
  % MeasureCoverage
    for c=1:24
      outfile = [OutDir{i} '/chr' num2str(c) '.cov'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
         banner = [short_sample 'MCV' tn(1) num2str(c)];
         cmd = ['"java -classpath ' java_classpath ' '...
               'MeasureCoverage ' BAMFile{i} ' ' BAIFile{i} ' ' outfile ' ' ...
                num2str(P.coverage_filter_flag) ' ' num2str(P.coverage_window_size) ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end % if P.MeasureCoverage

  if P.CountMQZ
    % CountMQZ
    for c=1:24
      outfile = [OutDir{i} '/chr' num2str(c) '.mqz'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'MQZ' tn(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' '...
                'CountMQZ ' BAMFile{i} ' ' BAIFile{i} ' ' outfile ' ' ...
                num2str(P.coverage_window_size) ' ' num2str(c) '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
  end % if P.CountMQZ

 end % next i   (T/N)

 if ~all_done
   fprintf('Waiting for dRangerPreprocess jobs to finish\n');
   bwait(jobs);
 end
end % while(~all_done)

if P.CatSortJoin

  % CONCATENATE

  jobs = [];
  all_done = false;
  while(~all_done)
    all_done = true;
    for i=1:2
      if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
      else if i==1, tn='sample';else break; end, end
      dbam = dir(BAMFile{i});
      outfile = [OutDir{i} '/all.weird'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'CAT' tn(1)];
        cmd = ['"cat ' OutDir{i} '/chr*.weird > ' outfile '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
    if ~all_done
      fprintf('Waiting for concatenate to finish\n');
      bwait(jobs);
    end
  end

  % SORT

  jobs = [];
  all_done = false;
  while(~all_done)
    all_done = true;
    for i=1:2
      if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
      else if i==1, tn='sample';else break; end, end
      dbam = dir(BAMFile{i});
      outfile = [OutDir{i} '/all.weird.sorted'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'SORT' tn(1)];
        cmd = ['"sort -T ' OutDir{i} ' ' OutDir{i} '/all.weird > ' outfile '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
    if ~all_done
      fprintf('Waiting for sort to finish\n');
      bwait(jobs);
    end
  end

  % JOIN

  jobs = [];
  all_done = false;
  while(~all_done)
    all_done = true;
    for i=1:2
      if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
      else if i==1, tn='sample';else break; end, end
      dbam = dir(BAMFile{i});
      outfile = [OutDir{i} '/all.weird.joined'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout)
        if dout.datenum>=dbam.datenum, uptodate = true;
        else fprintf('%s needs to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banner = [short_sample 'JOIN' tn(1)];
        cmd = ['"cat ' OutDir{i} '/all.weird.sorted | perl ' joinweird_perlscript ...
          ' > ' outfile '"'];
        jobs = [jobs;bsub(cmd,banner)];
      end
    end
    if ~all_done
      fprintf('Waiting for join to finish\n');
      bwait(jobs);
    end
  end
end % if P.CatSortJoin

if did_nothing
  fprintf('All dRanger input files already up-to-date\n');
else
  fprintf('dRangerPreprocess finished.\n');
end

catch me, excuse(me); end
