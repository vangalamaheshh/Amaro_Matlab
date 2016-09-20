function save2(X,fname)

if nargin~=2 || ischar(X) || ~ischar(fname), error('usage:  save2(X,fname)'); end
if ~isstruct(X), error('X should be a struct'); end

if exist(fname,'file') && ~exist(fname,'dir')
  error('%s already exists and is not a directory',fname);
end

ensure_dir_exists(fname);

f = fieldnames(X);
for i=1:length(f)
  tmp = getfield(X,f{i});
  save([fname '/' f{i} '.mat'],'tmp');
end

