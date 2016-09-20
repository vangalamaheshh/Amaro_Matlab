function bkill(jobid)

if ~iscell(jobid), jobid = {jobid}; end

fprintf('bkill the following jobs:\n');
for i=1:length(jobid)
  fprintf('%s\n',jobid{i});
end
a = input('Are you sure? (Y/N)','s');
if ~strcmpi(a,'Y')
  fprintf('Aborted bkill.\n');
  return
end

cmd = 'bsub';
for i=1:length(jobid)
  cmd = [cmd ' ' jobid{i}];
end
system(cmd);
 
