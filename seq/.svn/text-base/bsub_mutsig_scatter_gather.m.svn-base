function result = bsub_mutsig_scatter_gather(P)

if ~exist('P','var'), P=[]; end

if ~isfield(P,'covfile') && isfield(P,'terrfile'), P = rename_field(P,'terrfile','covfile'); end

P = impose_default_value(P,'outdir','*required*');
P = impose_default_value(P,'jobcount','*required*');
P = impose_default_value(P,'mutfile','*required*');
P = impose_default_value(P,'covfile','*required*');
P = impose_default_value(P,'build','*required*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = impose_default_value(P,'submit_jobs',true);
P = impose_default_value(P,'wait_for_jobs',true);
P = impose_default_value(P,'skip_gather',false);
P = impose_default_value(P,'assume_no_jobs_running_or_pending',false);
P = impose_default_value(P,'queue','cga');
P = impose_default_value(P,'bsub_queue',P.queue);
P = impose_default_value(P,'queue_gather',P.bsub_queue);
P = impose_default_value(P,'bsub_queue_gather',P.queue_gather);
P = impose_default_value(P,'queue_scatter',P.bsub_queue);
P = impose_default_value(P,'maxtries','3');
P = impose_default_value(P,'min_before_kill',inf);
if isnumeric(P.maxtries), P.maxtries = num2str(P.maxtries); end
P = impose_default_value(P,'bsub_queue_scatter',P.queue_scatter);
P = impose_default_value(P,'bsub_priority','60');
P = impose_default_value(P,'bsub_user_priority',P.bsub_priority);
if isnumeric(P.bsub_user_priority), P.bsub_user_priority = num2str(P.bsub_user_priority); end
P = impose_default_value(P,'bsub_project','Firehose');
P = impose_default_value(P,'bsub_mem','15');
if isnumeric(P.bsub_mem), P.bsub_mem = num2str(P.bsub_mem); end
P = impose_default_value(P,'check_matlab_version',true);
P = impose_default_value(P,'run_matlab_path','/cga/tcga-gsc/home/lawrence/bin/running_compiled_matlab/run_matlab.py');
P = impose_default_value(P,'pause_time_in_seconds',10);
P = impose_default_value(P,'matlab_use_display',true);
P = impose_default_value(P,'java_opts','-Xmx2g');
P = impose_default_value(P,'bsub_exclude','"EXCLUDE(167)"');
P = impose_default_value(P,'banner','MutSig');
P = impose_default_value(P,'mutsig_version','latest');
P = impose_default_value(P,'catfile',[]);
if isfield(P,'patlist') && ~isfield(P,'patfile'), P = rename_field(P,'patlist','patfile'); end
P = impose_default_value(P,'patfile',[]);
if isfield(P,'genelist') && ~isfield(P,'genefile'), P = rename_field(P,'genelist','genefile'); end
P = impose_default_value(P,'genefile',[]); % /cga/tcga-gsc/home/lawrence/capture/Refseq_exons_good_20101221_genelist.txt
P = impose_default_value(P,'individual_set_id',[]);
%P = impose_default_value(P,'mutation_preprocessing_report_file',[]);
P = impose_default_value(P,'cosmic_file','/xchip/cga1/annotation/db/cosmic/CosmicMutantExport_v48_260710_p25.tsv');
P = impose_default_value(P,'geneset_collection_file','/xchip/cga1/annotation/db/genesets/gsea_canonical_pathway_plus_HMT67.txt');
if strcmp(P.build,'hg18')
  P = impose_default_value(P,'refseq_file','/xchip/cga1/annotation/db/ucsc/hg18/R.mat');
  P = impose_default_value(P,'build_dir','/xchip/cga1/annotation/db/ucsc/hg18');
elseif strcmp(P.build,'hg19')
  P = impose_default_value(P,'refseq_file','/xchip/cga1/annotation/db/ucsc/hg19/R.mat');
  P = impose_default_value(P,'build_dir','/xchip/cga1/annotation/db/ucsc/hg19');
else
  error('unknown build %s',P.build);
end

if isfield(P,'params_file')
  error('Please do not specify P.params_file\n\t--this is handled by bsub_mutsig_scatter_gather');
end

% make sure there are no missing input files
flds = grepi('file',fieldnames(P));
if ~isempty(flds)
%  fprintf('Confirming the existence of the following files:\n'); disp(flds);
  for i=1:length(flds), if ~isempty(P.(flds{i})), demand_file(getfield(P,flds{i})), end, end
end
if isfield(P,'summed_cov_track') && ~isempty(P.summed_cov_track)
  demand_file(P.summed_cov_track);
  demand_file(regexprep(P.summed_cov_track,'\.fwb$','.fwi'));
end

% check jobcount
if ~isnumeric(P.jobcount) || P.jobcount<1 || isinf(P.jobcount) || isnan(P.jobcount)
  error('invalid P.jobcount');
end
% make sure jobcount<=length(genelist)
if isfield(P,'genefile') && exist(P.genefile,'file')
  tmp = load_lines(P.genefile);
  tmp = grepv('^gene|name',tmp);
  if P.jobcount>length(tmp)
    fprintf('Decreasing jobcount to length(genelist) = %d\n',length(tmp));
    P.jobcount = length(tmp);
  end
end
if P.jobcount~=round(P.jobcount)
  fprintf('Rounding P.jobcount up to %d\n',ceil(P.jobcount));
  P.jobcount = ceil(P.jobcount);
end

%% make sure we're using right version of matlab
if P.check_matlab_version
  [result output] = system('which matlab');
  if (result~=0), error('Error checking which matlab'); end
  if contains(output,'matlab_2009')
    % OK
  elseif contains(output,'matlab_2010')
    fprintf('Currently using matlab 2010\n');
    fprintf('Need to switch to matlab 2009.\n');
    fprintf('Please exit Matlab and do\n');
    fprintf('       reuse Matlab-2009a\n');
    fprintf('Then try again.\n');
    error('Wrong Matlab version');
  else
    fprintf('Not sure if this will work with the version of Matlab in use, which is:  ');
    disp(output);
  end
end

% find directory with latest mutsig version
libdir = '';
if ischar(P.mutsig_version)
  if strcmpi(P.mutsig_version,'latest')
    P.mutsig_version = inf;
  elseif P.mutsig_version(1)=='/'
    libdir = P.mutsig_version;
  else
    tmp = str2double(P.mutsig_version);
    if ~isnan(tmp)
      P.mutsig_version = tmp; 
    else
      error('invalid P.mutsig_version');
    end
  end
end
if isnumeric(P.mutsig_version)
  tl = '/xchip/tcga/gdac_prod/applications/genepattern/gp-3.2.4-9272/Tomcat/../../taskLib';
  v=[]; v.dir = direc([tl '/MutSigRun.*']);
  v = parse_in(v,'dir','MutSigRun.(\d+).(\d+)$',{'ver','svn'},1:2);
  if isempty(v), error('Could not find MutSigRun in task directory'); end
  if isinf(P.mutsig_version)
    [tmp vidx] = max(v.ver);
    if isempty(vidx), error('Error locating latest version of MutSigRun'); end
  else
    vidx = find(v.ver==P.mutsig_version);
    if isempty(vidx)
      fprintf('Could not find MutSigRun version:  '); disp(P.mutsig_version);
      error('Leave P.mutsig_version unspecified to choose latest version.');
    end
  end
  libdir = v.dir{vidx};
elseif ischar(P.mutsig_version)
  v = parse(P.mutsig_version,'MutSigRun.(\d+).(\d+)$',{'ver','svn'},1:2);
  if ~isempty(v.ver)
    vidx = v.ver;
  end
else
  error('unknown format for P.mutsig_version');
end
if isempty(libdir)
  error('unvalid P.mutsig_version');
end
libexe = [libdir '/fh_MutSigRun'];
demand_file(libexe);
[tmp attrib] = fileattrib(libexe);
if ~attrib.UserExecute
  if exist('vidx','var') && length(v.ver)>=vidx
    error('MutSigRun v%d needs a chmod +x from Firehose.',v.ver(vidx));
  else
    error('%s needs a chmod +x from Firehose.',libexe);
  end
end

% create output directory
ensure_dir_exists(P.outdir);

% write parameters file
params_file = [P.outdir '/params.txt'];
write_params_file(P,params_file);

% BSUB SCRIPT

basescript = sprintf([...
    '#BSUB -r\n'...
    '#BSUB -R "rusage[mem=' P.bsub_mem ']"\n'...
    '#BSUB -R "select[tmp>1000 && scratch>1]"\n'...
    '#BSUB -q <QUEUE>\n'...
    '#BSUB -P ' P.bsub_project '\n'...
    '#BSUB -sp ' P.bsub_user_priority '\n'...
    '#BSUB -Q ' P.bsub_exclude '\n'...
    '#BSUB -cwd <JOBDIR>\n'...
    '#BSUB -e job.%%J.err\n'...
    '#BSUB -o job.%%J.out\n'...
    '#BSUB -J ' P.banner '<WHICHJOB>\n'...
    '#BSUB -E "rm -f status_pending && rm -f *.out && rm -f *.err && echo ' P.java_opts ' > java.opts && touch status_running"\n'...
    'tries=`cat <TRIESFILE>`\n'...
    'echo -n `expr ${tries} + 1` > <TRIESFILE>\n'...
    'python ' P.run_matlab_path ' '...
]);

if P.matlab_use_display, basescript = [basescript ' --with-display ']; end

basescript = [basescript ...
       libdir ' fh_MutSigRun '];

%if vidx>151      % "libdir" was added as first parameter sometime after v151 and before v168
  basescript = [basescript ' ' libdir ' '];
%end

basescript = [basescript ...
       ' -jobcount ' num2str(P.jobcount) ' -whichjob <WHICHJOB> ' ...
       ' -b ' P.build ' -bd ' P.build_dir ' -maf ' P.mutfile ' -cov ' P.covfile ...
       ' -cosmic ' P.cosmic_file ' -geneset ' P.geneset_collection_file ' -refseq ' P.refseq_file ' ' ... 
       ' -p ' params_file ' ' ...
];

opt1 = {'individual_set_id','catfile','patfile','genefile'};
opt2 = {'-i'               ,'-cat'   ,'-pat'   ,'-genes'};
for i=1:length(opt1)
  tmp = getfield(P,opt1{i});
  if ~isempty(tmp)
    basescript = [basescript ' ' opt2{i} ' ' tmp ' '];
  end
end

basescript = [basescript sprintf([...
    '\n'...
    'return_code=$?\n'...
    'if [ -e core.* ]; then\n'...
      'return_code=999\n'...
    'fi\n'...
    'if [ ${return_code} -ne 0 ]; then\n'...
      'tries=`cat <TRIESFILE>`\n'...
      'if [ ${return_code} -ne 130 -a ${tries} -lt <MAXTRIES> ]; then\n'...
        'rm status_running\n'...
        'exit 167;\n'...
      'else\n'...
        'echo ${return_code} > status_exited;\n'...
        'rm status_running\n'...
        'exit 1;\n'...
       'fi\n'...
    'else\n'...
      'touch status_success\n'...
      'rm status_running\n'...
    'fi\n'...
])];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTER

firstpass = true;
finished = false;
while(~finished)

  if P.assume_no_jobs_running_or_pending
    fprintf('Assuming no jobs running/pending\n');
  end

  nsuccess=0;
  nexited=0;
  npending=0;
  nrunning=0;
  nsubmit=0;
  batchsubmit = cell(P.jobcount,1);
  for whichjob=1:P.jobcount
%    disp(whichjob)
    jobname = ['scatter.' sprintf('%010.0f',whichjob)];
    jobdir = [P.outdir '/' jobname];
    need_to_submit = true;
    if exist(jobdir,'dir')
      for tryno=1:(3-2*(firstpass==true))
        if exist([jobdir '/status_success'],'file'), nsuccess=nsuccess+1; need_to_submit=false; break; end
        if exist([jobdir '/status_exited'],'file'), nexited=nexited+1; need_to_submit = false; break; end
        d = direc([jobdir '/core.*']);
        if ~isempty(d) % job apparently crashed with core dump
          fprintf('Job %d core dump\n',whichjob);
          delete([jobdir '/core.*']); % delete core dump file (these add up to a lot of disk space)
          if exist([jobdir '/status_running'],'file'), delete([jobdir '/status_running']); end
          save_textfile('core',[jobdir '/status_exited']); nexited=nexited+1; need_to_submit = false; break
        end
        if ~P.assume_no_jobs_running_or_pending
          if exist([jobdir '/status_running'],'file')
            d = dir([jobdir '/status_running']);
            mins_running = (now-d.datenum)*24*60;
            if mins_running >= P.min_before_kill
              % job has been "running" too long: time to kill+restart
              % (for now, will just resubmit)
              fprintf('Job %d has been running too long: resubmitting.\n',whichjob);
            else
              nrunning=nrunning+1; need_to_submit = false; break;
            end
          end
          if exist([jobdir '/status_pending'],'file'), npending=npending+1; need_to_submit = false; break; end
        end
        % else: no indication what happened to this job
        if (firstpass==false) pause(2); end % pause and try checking again (to avoid problem of delayed appearance of "status_success" file)
      end % next tryno
      if need_to_submit && ~isempty(direc([jobdir '/*.out']))  % job was apparently killed after producing some output: call it "exited"
        fprintf('Assuming job %d exited\n',whichjob);
        save_textfile('killed',[jobdir '/status_exited']); nexited=nexited+1;
        need_to_submit = false;
      end
    end
    if need_to_submit
      nsubmit=nsubmit+1;
      ensure_dir_exists(jobdir);
      triesfile = [P.outdir '/' jobname '.tries'];
      save_textfile('0',triesfile);
      scriptfile = [P.outdir '/' jobname '.script'];
      script = basescript;
      script = regexprep(script,'<QUEUE>',P.bsub_queue_scatter);
      script = regexprep(script,'<JOBDIR>',jobdir);
      script = regexprep(script,'<WHICHJOB>',num2str(whichjob));
      script = regexprep(script,'<TRIESFILE>',triesfile);
      script = regexprep(script,'<MAXTRIES>',P.maxtries);
      save_textfile(script,scriptfile);
      batchsubmit{whichjob} = sprintf('bsub < %s && touch %s/status_pending\n', scriptfile, jobdir);
    end
  end
  
  status = '[scatter]   ';
  if nsubmit>0,  status = [status sprintf('%d SUBMIT   ',nsubmit)]; end
  if npending>0, status = [status sprintf('%d PEND     ',npending)]; end
  if nrunning>0, status = [status sprintf('%d RUN      ',nrunning)]; end
  if nexited>0,  status = [status sprintf('%d EXITED   ',nexited)]; end
  if nsuccess>0, status = [status sprintf('%d SUCCESS  ',nsuccess)]; end
  if exist('oldstatus','var')
    if strcmp(oldstatus,status), fprintf('\r'); else fprintf('\n'); end
  end
  fprintf('%s   %s', datestr(now),status);
  oldstatus = status;

  if nsubmit>0
    batchsubmit = cat(2,batchsubmit{:});
    batchsubmitfile = [P.outdir '/scatter.submit.script'];
    save_textfile(batchsubmit,batchsubmitfile);
    if P.submit_jobs
      [result output] = unix(['source ' batchsubmitfile]);
      if result~=0
        fprintf('\nLSF: nonzero return code on batchsubmit\n');
        disp(output);
      end
    end
  end

  finished = (nexited+nsuccess == P.jobcount);

  if ~finished && ~P.wait_for_jobs, fprintf('\n'); result=2; return; end

  if (~finished) pause(P.pause_time_in_seconds); end
  firstpass=false;
end

if nexited>0
  fprintf('\nSome scatter jobs failed.\n');
  result=-1;
  return
end

if P.skip_gather
  result=0;
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GATHER

finished = false;
while(~finished)

  if P.assume_no_jobs_running_or_pending
    fprintf('Assuming no jobs running/pending\n');
  end

  nsuccess=0;
  nexited=0;
  npending=0;
  nrunning=0;
  nsubmit=0;

  jobname = 'gather';
  jobdir = P.outdir;
  need_to_submit = true;
  if exist(jobdir,'dir')
    for tryno=1:3
      if exist([jobdir '/status_success'],'file'), nsuccess=nsuccess+1; need_to_submit=false; break; end
      if exist([jobdir '/status_exited'],'file'), nexited=nexited+1; need_to_submit = false; break; end
      if ~isempty(direc([jobdir '/core.*'])) % job apparently crashed with core dump
        save_textfile('core',[jobdir '/status_exited']); nexited=nexited+1; need_to_submit = false; break
      end
      if ~P.assume_no_jobs_running_or_pending
        if exist([jobdir '/status_running'],'file'), nrunning=nrunning+1; need_to_submit = false; break; end
        if exist([jobdir '/status_pending'],'file'), npending=npending+1; need_to_submit = false; break; end
      end
      % else: no indication what happened to this job
      pause(2);  % pause and try checking again (to avoid problem of delayed appearance of "status_success" file)
    end % next tryno
    if need_to_submit && ~isempty(direc([jobdir '/*.out']))  % job was apparently killed after producing some output: call it "exited"
      save_textfile('killed',[jobdir '/status_exited']); nexited=nexited+1;
      need_to_submit = false;
    end
  end
  if need_to_submit
    nsubmit=nsubmit+1;
    triesfile = [P.outdir '/' jobname '.tries'];
    save_textfile('0',triesfile);
    scriptfile = [P.outdir '/' jobname '.script'];
    script = basescript;
    script = regexprep(script,'<QUEUE>',P.bsub_queue_scatter);
    script = regexprep(script,'<JOBDIR>',jobdir);
    script = regexprep(script,'<WHICHJOB>','0');
    script = regexprep(script,'<TRIESFILE>',triesfile);
    script = regexprep(script,'<MAXTRIES>','1');   % do not increase for gather: will delete files in P.outdir!
    save_textfile(script,scriptfile);
    if P.submit_jobs
      submit = sprintf('bsub < %s && touch %s/status_pending\n', scriptfile, jobdir);
    end
  end
  
  status = '[gather]    ';
  if nsubmit>0,  status = [status sprintf('%d SUBMIT   ',nsubmit)]; end
  if npending>0, status = [status sprintf('%d PEND     ',npending)]; end
  if nrunning>0, status = [status sprintf('%d RUN      ', nrunning)]; end
  if nexited>0,  status = [status sprintf('%d EXITED   ',nexited)]; end
  if nsuccess>0, status = [status sprintf('%d SUCCESS  ',nsuccess)]; end
  if exist('oldstatus','var')
    if strcmp(oldstatus,status), fprintf('\r'); else fprintf('\n'); end
  end
  fprintf('%s   %s', datestr(now),status);
  oldstatus = status;

  if nsubmit>0
    submitfile = [P.outdir '/gather.submit.script'];
    save_textfile(submit,submitfile);
    [result output] = unix(['source ' submitfile]);
    if result~=0
      fprintf('\nLSF: nonzero return code on submit\n');
      disp(output);
    end
  end

  finished = (nexited+nsuccess == 1);

  if ~finished && ~P.wait_for_jobs, fprintf('\n'); result=3; return; end

  if (~finished) pause(P.pause_time_in_seconds); end
end

if nexited>0
  fprintf('\nGather job failed.\n');
else
  fprintf('\nGather job finished successfully.\n');
end

result=1;
