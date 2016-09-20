function segseq_preprocess(sample,P)
% segseq_preprocess(sample,P)

% Mike Lawrence 2009, updated 2010 for absolute-path mode

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'unpaired_flag',false);
P=impose_default_value(P,'blacklist','none');
P=impose_default_value(P,'bsub_queue','cga');
P=impose_default_value(P,'queue',P.bsub_queue);

if isempty(sample)    % absolute-path mode (not using "lawrence" directory)
  fprintf('Using absolute-path mode\n');
  P=impose_default_value(P,'tumor_bam','*required*');
  P=impose_default_value(P,'normal_bam','*required*');
  P=impose_default_value(P,'tumor_output_dir','*required*');
  P=impose_default_value(P,'normal_output_dir','*required*');
  P=impose_default_value(P,'sample_name','');
  if ~iscell(P.tumor_bam)
    nsamps=1;
    if iscell(P.normal_bam) || iscell(P.tumor_output_dir) || iscell(P.normal_output_dir) || iscell(P.sample_name)
      error('All or none of the inputs must be cell arrays');
    end
    P.BAMFile = {P.tumor_bam P.normal_bam};
    P.OutDir = {P.tumor_output_dir P.normal_output_dir};
    short_sample = {P.sample_name};
  else
    nsamps = length(P.tumor_bam);
    if length(P.normal_bam)~=nsamps || length(P.tumor_output_dir)~=nsamps || length(P.normal_output_dir)~=nsamps || ...
          length(P.sample_name)~=nsamps, error('All inputs must be same length');
    end
    BAMFile = [P.tumor_bam P.normal_bam];
    OutDir = [P.tumor_output_dir P.normal_output_dir];
    ensure_dir_exists(OutDir);
    short_sample = P.sample_name;
  end
else   % "classic" mode
  fprintf('Using classic file-addressing mode\n');
  if iscell(sample), error('Multiple samples not supported in "classic" mode'); end
  for i=1:2
    if i==1, tn='normal';else tn='tumor'; end
    BAMFile{i,1} = ['/xchip/cga1/lawrence/' sample '/' tn '.bam'];
    OutDir{i,1} = ['/xchip/cga1/lawrence/' sample '/' tn '_SS'];
    ensure_dir_exists(OutDir{i});
  end
  short_sample = sample_to_short_sample(sample);
end

try

java_classpath = get_jcp;

SSP_class = 'org/broadinstitute/cga/tools/seq/SegSeqPreprocess';
MLL_class = 'org/broadinstitute/cga/tools/seq/MakeLanelist';

did_nothing = true;
all_done = false;
while(~all_done)
  all_done = true;
  jobs=[];
  cmds = {}; banners = {};
  for sampno=1:nsamps
   for i=1:2
    if i==1, tn='normal';else tn='tumor'; end
    dbam = dir(BAMFile{sampno,i});
    if isempty(dbam), error('%s not found',BAMFile{sampno,i}); end
    for c=1:24
      outfile = [OutDir{sampno,i} '/chr' num2str(c) '.txt'];
      dout = dir(outfile);
      uptodate = false;
      if ~isempty(dout) & dout.bytes>1000
        if dout.datenum >= dbam.datenum, uptodate = true;
        else fprintf('%s need to be refreshed:\n',outfile); end
      end
      if ~uptodate
        did_nothing = false;
        all_done = false;
        banners{end+1,1} = [short_sample{sampno} 'SSP' tn(1) num2str(c)];
        cmd = ['"java -classpath ' java_classpath ' ' SSP_class ' ' ...
           BAMFile{sampno,i} ' ' P.blacklist ' ' outfile ' ' num2str(c) ' mqual'];
        if P.unpaired_flag, cmd = [cmd ' unpaired']; end
        cmd = [cmd '"'];
        cmds{end+1,1} = cmd;
        % submit jobs as they accumulate
        if length(cmds)>=60
          jobs = [jobs; bsub(cmds,banners,P)];
          cmds = {}; banners = {};
        end
      end
    end
    outfile = [OutDir{sampno,i} '/lanelist.txt'];
    dout = dir(outfile);
    if ~isempty(dout) & dout.bytes>0
      if dout.datenum >= dbam.datenum, uptodate = true;
      else fprintf('%s need to be refreshed:\n',outfile); end
    end
    if ~uptodate
      did_nothing = false;
      all_done = false;
      banners{end+1,1} = [short_sample{sampno} 'LANES' tn(1)];
      cmds{end+1,1} = ['"java -classpath ' java_classpath ' ' MLL_class ' ' ...
         BAMFile{sampno,i} ' ' outfile '"'];
      % submit jobs as they accumulate
      if length(cmds)>=60
        jobs = [jobs; bsub(cmds,banners,P)];
        cmds = {}; banners = {};
      end
    end
   end % next i   (T/N)
  end % next sampno

  if ~isempty(cmds)
    jobs = [jobs; bsub(cmds,banners,P)];
    cmds = {}; banners = {};
  end

  if ~all_done
    fprintf('Waiting for SegSeqPreprocess to finish\n');
    bwait(jobs);
  end
end % while(~all_done)

if did_nothing
  fprintf('All SegSeq input files already up-to-date.\n');
else
  fprintf('All files now exist.\n');
end

catch me, excuse(me); end
