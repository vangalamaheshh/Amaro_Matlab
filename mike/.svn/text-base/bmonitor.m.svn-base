function bmonitor(secs)
% Mike Lawrence 2009-006-30

if ~exist('secs','var'), secs = 60; end

while 1
    j = bjobs;
    if slength(j)>0
      j.type = regexprep(j.job_name,'[a-z0-9_]','');
      xcount(j.type,j.stat);
    else
      fprintf('No jobs.\n');
    end
    pause(secs);
end
