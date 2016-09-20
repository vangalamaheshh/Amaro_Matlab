function globalrep(basedir,recursive,files,pattern,replacement,prompt)
%
% globalrep(files,pattern,replacement)
%
% basedir = directory to start from
% recursive = true/false, navigate downward through directories
% files = files to do search-and-replace on
% pattern,replacement = for passing to regexprep
% prompt = true/false, prompt on each replacement (default=true)
%

if ~exist('prompt','var'), prompt=true; end

if iscell(basedir)
  for i=1:length(basedir)
    globalrep(basedir{i},recursive,files,pattern,replacement,prompt);
  end
  return
end

if ~isempty(basedir) && basedir(end)=='/', basedir(end)=[]; end  % trim trailing forwardslash

fprintf('Checking directory %s\n',basedir);

% first handle files
d = dir([basedir '/' files]);
for i=1:length(d)
  if d(i).isdir, continue; end
  if d(i).bytes > 300000000, fprintf('Skipping %s, over 300MB\n',d(i).name), continue; end
  fname = [basedir '/' d(i).name];
  fprintf('Checking file %s\n',fname);
  X = load_textfile(fname);
  any_reps = false;
  startpos = 1;
  while(startpos<length(X))
    S = regexp(X(startpos:end),pattern,'start')+startpos-1;
    if isempty(S), break; end
    if startpos==1
      fprintf('\n==========================================================\n');
      fprintf('%d matches found in file:\n%s\n',length(S),fname);
      fprintf('==========================================================\n');
    end
    E = regexp(X(startpos:end),pattern,'end')+startpos-1;
    tmp = find(X==char(10))';
    L = [[1;tmp(1:end-1)+1] tmp];
    l = find(S(1)>=L(:,1) & S(1)<=L(:,2));
    ln = X(L(l,1):L(l,2));
    fprintf('\n%s',ln);
    r = regexprep(X(S(1):E(1)),pattern,replacement);
    ln = [ln(1:S(1)-L(l,1)) r ln(E(1)-L(l,1)+2:end)]; 
    fprintf('===[to]===>\n');
    fprintf('%s',ln);
    if prompt
      answ = [];
      while(~strcmpi(answ,'y')&&~strcmpi(answ,'n'));
        answ = input('[REPLACE? (y/n)] ','s');
      end
    else
      answ = 'y';
    end
    if strcmpi(answ,'y')
      X = [X(1:(S(1)-1)) r X(E(1)+1:end)];
      any_reps = true;
      startpos = S(1) + length(r);
    else
      startpos = E(1) + 1;
    end
  end
  if any_reps, save_textfile(X,fname); end
end

% then handle directories
if recursive
  d = dir(basedir);
  for i=1:length(d)
    if ~d(i).isdir, continue; end
    if strcmp(d(i).name,'.') || strcmp(d(i).name,'..'), continue; end

    % DIRECTORIES TO AVOID (temporary measure)
    if strcmp(d(i).name,'bak'), continue; end     % TEMPORARY
    if strncmp(d(i).name,'maps',4), continue; end % TEMPORARY
    if strncmp(d(i).name,'allfigs',7), continue; end % TEMPORARY
    if strncmp(d(i).name,'newfigs',7), continue; end % TEMPORARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    globalrep([basedir '/' d(i).name],recursive,files,pattern,replacement,prompt);
  end
end
