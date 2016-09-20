function db = load_dbSNP(build)
% load_dbSNP(build)
%
% loads the specified dbSNP build and returns it as a cell(24,1) array of "pos" arrays
%   for dbSNP positions, it's the rs______ number of the dbSNP
%   for non-dbSNP positions, it's 0.
%
% Mike Lawrence 2011-09-29

if ~exist('build','var')
  fprintf('Assuming build hg18\n');
  build='hg18';
end

if isnumeric(build)
  if build==18 || build==130
    build = '130';
  elseif build==19 || build==132
    build = '132';
  else
    error('unknown build %d',build);
  end
end

if ~ischar(build), error('build should be string'); end

if strcmp(build,'hg18') || strcmp(build,'130')
  dbSNPdir = '/xchip/cga1/lawrence/db/hg18/dbsnp/130';
elseif strcmp(build,'hg19') || strcmp(build,'132')
  dbSNPdir = '/xchip/cga1/lawrence/db/hg19/dbsnp/132';
else
  error('unknown build:  ');disp(build);
end

db = cell(24,1);
for chr=1:24
  fprintf('chr%d ',chr);
  fname = [dbSNPdir '/chr' num2str(chr) '.dbSNP'];
  d = dir(fname);
  if isempty(d), error('Not found: %s',fname); end
  len = floor(d.bytes/4);
  db{chr} = get_block(fname,'int',0,len-1);
end, fprintf('\n');

