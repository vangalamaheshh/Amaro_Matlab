function bwait(jobs,n_confirms_needed)
% Mike Lawrence 2009

if isempty(jobs), return; end

if ischar(jobs), jobs={jobs}; end

if ~exist('n_confirms_needed','var'), n_confirms_needed=1; end

nwaitingsave = inf; ujsave = []; usave = [];

max_problems = 3;
problems = 0;
confirms = 0;
while 1
  try
    nwaiting = nan;
    nrunning = nan;
    j = bjobs;
    nrunning = length(j.jobid);
    if nrunning>0
      [jobs_left jidx] = intersect(j.jobid,jobs);
      nwaiting = length(jobs_left);
    else
      nwaiting = 0;
    end
  catch me
    fprintf('Problem with bwait:\n');
    problems = problems + 1;
    if problems>max_problems, keyboard; end
    me
  end
  if nwaiting>0
    [u ui uj] = unique(j.stat(jidx));
  else
    u = [];
    ui = [];
    uj = [];
  end
  if nwaiting~=nwaitingsave | length(uj)~=length(ujsave) | any(uj~=ujsave) |...
             length(usave)~=length(u) | any(~strcmp(usave,u))
    nwaitingsave = nwaiting;
    ujsave = uj;
    usave = u;
    fprintf('%s   Waiting for %d jobs',datestr(now),nwaitingsave);
    if nwaiting>0
      fprintf(':');
      [u ui uj] = unique(j.stat(jidx));
      for i=1:length(u)
        fprintf(' %s(%d)',u{i},sum(uj==i));
      end
    end
    fprintf('\n');
  end
  if nwaiting>0 && nwaiting<4
    fprintf('\nHere is the output of the last %d job(s) we''re waiting for:\n',nwaiting);
    for i=1:length(jobs_left)
      [status result] = system(['bpeek ' jobs_left{i}]);
      fprintf('  JOB %s:\n',jobs_left{i});
      result = regexprep(result,'<< output from stdout >>','');
      lines = text_to_lines(result);
      nl = length(lines);
      num_to_show = min(3,nl);
      for j=num_to_show:-1:1
        fprintf('    %s\n',lines{end-j+1});
      end
    end
  end
  if nwaiting==0
    confirms=confirms+1;
    if confirms==n_confirms_needed
      break
    else
      pause(1); 
      continue
    end
  else
    confirms = 0;
    pause(15);
  end
end
