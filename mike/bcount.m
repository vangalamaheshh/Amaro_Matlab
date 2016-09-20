function nwaiting = bcount(jobs)
% counts how many of the jobs are not finished
% Mike Lawrence 2009-09-25

if isempty(jobs), return; end

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
    fprintf('Problem with bcount:\n');
    me
  end
  if ~isnan(nwaiting), break; end
end
