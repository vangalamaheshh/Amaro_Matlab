function X = load2(fname)

if nargin~=1 || nargout~=1 || ~ischar(fname), error('usage:  X = load2(fname)'); end

if ~exist(fname,'dir')
  error('directory %s does not exist',fname);
end

X = [];
d.file = direc([fname '/*.mat']);
d = parse_in(d,'file',[fname '/(.*)\.mat$'],'field');
for i=1:slength(d)
  clear tmp;
  load(d.file{i},'tmp');
  if ~exist('tmp','var'), error('%s does not contain "tmp" variable',d.file{i}); end
  X = setfield(X,d.field{i},tmp);
end
