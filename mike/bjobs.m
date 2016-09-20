function J = bjobs()

while(1)

try
  [result output] = unix('bjobs');

  output = split(output,char(10));

  keep=[];
  for i=1:length(output)
    if isempty(output{i}), continue; end
    if contains(output{i},'Cannot connect'), continue; end
    if contains(output{i},'Please wait'), continue; end
    keep=[keep;i];
  end
  output = output(keep);

  if ~isempty(grep('please',output,1))
    fprintf('bjobs is still not working properly\n');
   keyboard;
  end 

  if strcmp(output{1},'No unfinished job found')
    J = struct('jobid',[]);
    return;
  end

  max_col = max(cellfun('length',output));

  header = lower(output{1});
  jobs = output(2:end);

  start_col = find([32 header]==32 & [header 0]~=32);
  end_col = [start_col(2:end)-1 max_col];
  nc = length(start_col);

  fields = cell(nc,1);
  for c=1:nc, fields{c} = remove_trailing_whitespace(header(start_col(c):min(end_col(c),length(header)))); end
  fields = genfieldname(fields);

  J = [];
  for c=1:nc
    x = cell(length(jobs),1);
    for i=1:length(jobs),x{i} = jobs{i}(start_col(c):min(end_col(c),length(jobs{i}))); end
    J = setfield(J,fields{c},remove_trailing_whitespace(x));
  end

  zzz = fieldnames(J);

  if ~isfield(J,'jobid')
    fn = fieldnames(J);
    idx = find(grepm('^.*jobid.*$',fn));
    % filter out stray keystrokes during the system call
    if ~isempty(idx)
      idx = idx(1);
      J = rename_fields(J,fn{idx},'jobid');
    end
  end

  if ~isfield(J,'jobid')
    fprintf('Problem with bjobs()\n');
    pause(5);
    J
  else
    return
  end

catch me
  fprintf('Poblem with bjobs()\n');
  pause(5);
  me
end

end
