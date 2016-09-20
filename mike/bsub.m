function job_no = bsub(cmds,banners,P)
% Mike Lawrence 2008-2012

if isempty(cmds)
  job_no = [];
  fprintf('No commands\n');
  return
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'logfile_outdir','/cga/tcga-gsc/home/lawrence/lsf');

try
  ede(P.logfile_outdir);
  tmpfile = [P.logfile_outdir '/lsfout.test.txt'];
  save_textfile('test',tmpfile);
%  delete(tmpfile);   % causes matlab to hang if there are too many files in the directory
catch me
  fprintf('Warning: Could not write to logfile_outdir %s\n',P.logfile_outdir);
  fprintf('         Will use write instead to /tmp/lsf\n');
  P.logfile_outdir = '/tmp/lsf';
  try
    ede(P.logfile_outdir);
  catch me
    error('Failed to create /tmp/lsf');
  end
end 

if isfield(P,'bsub_queue'), error('name of field is "queue" not "bsub_queue"'); end

P = impose_default_value(P,'queue','hour');

if ~iscell(cmds)
  cmds = {cmds};
  if exist('banners','var'), banners = {banners}; end
end

if ~exist('banners','var') || isempty(banners)
  banners = str2cell(sprintf('job%d\n',1:length(cmds)));
end

if strncmp(P.queue,'manual',6)
  if length(P.queue)<8, error('for manual queueing, specify "manual:/path/file.txt"'); end
  outfile = P.queue(8:end);
  try
    save_lines(cmds,outfile);
  catch me
    error('Error writing list of commandlines to %s',outfile);
  end
  fprintf('Saved list of job commandlines to %s\n\n',outfile);
  fprintf('Please submit these jobs to LSF and wait for them to finish.\n');
  input('When all jobs have finished, press ENTER to continue.\n');
  fprintf('Thank you so much!\n');
  job_no = [];
  return;
end

calls = [];
%calls = [calls 'reuse -q Java-1.6;'];   % needed since transition to CentOS  (but should be omitted with Solexa queue)
for i=1:length(cmds)
  call = 'bsub ';
  call = [call '-E "sleep 10;cd /cga/tcga-gsc;cd /seq/picard_aggregation" '];
  call = [call '-P cancer '];
  call = [call '-q ' P.queue ' '];
  call = [call '-o ' P.logfile_outdir '/%J.out '];
  call = [call '-r '];   % re-runnable (in case job gets SSUSP'ed)
%  call = [call '-R "rusage[thumper12_io=10]" '];
  call = [call '-J ' banners{i} ' '];
  call = [call cmds{i}];
  fprintf('\nSubmitting job: %s\n',call);
  if isempty(calls), calls = call; else calls = [calls ';' call]; end
end


[result output] = unix(calls);
if result~=0
  fprintf('bsub failed with exit code %d\n',result);
  fprintf('\nbsub output:\n\n');
  disp(output);  
  error('bsub failure');
end

output = split(output,char(10));
tmp = parse(output,'Job <(\d+)> is submitted',{'j'});
job_no = tmp.j(~cellfun('isempty',tmp.j));

if length(job_no) ~= length(cmds)
  error('Submitted %d jobs but only %d were accepted\n', length(cmds), length(job_no));
end
