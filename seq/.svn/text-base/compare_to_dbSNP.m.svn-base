function rs = compare_to_dbSNP(chr,pos,build)
% compare_to_dbSNP(chr,pos,build)
%
% given "chr" and "pos", which are a list of chr:pos genomic coordinates,
% returns a list of the same length:
%   for dbSNP positions, it's the rs______ number of the dbSNP
%   for non-dbSNP positions, it's 0.
%
% Mike Lawrence 2009-10-08

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

if ~isnumeric(chr), chr = convert_chr(chr); end
if length(chr)==1 & length(pos)>1, chr = repmat(chr,size(pos)); end
if length(chr)~=length(pos), error('chr and pos must be same length'); end

[c ci cj] = unique(chr);
if any(c<1 | c>24 | isnan(c)) error('chr should be 1-24'); end

rs = zeros(length(chr),1);

for i=1:length(c)
  chr = c(i); fprintf('chr%d ',chr);
  idx = find(cj==i);
  fname = [dbSNPdir '/chr' num2str(chr) '.dbSNP'];
  d = dir(fname);
  if isempty(d), error('Not found: %s',fname); end
  len = floor(d.bytes/4);
  idx = idx(pos(idx)<=len);
  if ~isempty(idx)
    rs(idx) = get_block(fname,'int',pos(idx)-1);   % FIXED BUG 2011-07-15:  added -1
  end
end, fprintf('\n');


