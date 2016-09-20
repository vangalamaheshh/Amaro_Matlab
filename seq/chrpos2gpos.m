function gpos = chrpos2gpos(chr,pos,P)
% gpos = chrpos2gpos(chr,pos)
%
% given a list of chromosomes <chr> and positions on them <pos>,
% computes <gpos>, the "genomic position", from 1 to ~3.2e9.
%
% uses hard-coded list of chromosome lengths
%
% Mike Lawrence 2009-03-17


if ~exist('P','var'), P=[]; end

if ischar(P)
  build=P;
  P=[];
end

if ~exist('build','var')
  build = [];
end

P = impose_default_value(P,'build',build);

if ~isempty(P.build)
  clen = load_chrlen(P.build) * 1.1;
else
  fprintf('Assuming human\n');
  clen = [...
      250000000   250000000   200000000   200000000   190000000   180000000   160000000   150000000   150000000 ...
      140000000   140000000   140000000   120000000   110000000   120000000   100000000    90000000    80000000 ...
      70000000    70000000    60000000    60000000   160000000    60000000 ...
         ]';
end

ct = length(clen);

% error checking on chr

if ~isnumeric(chr) error('chr must be numeric'); end
if size(chr,1)==1, chr=chr'; end
if any(chr-round(chr)) error('chr must be whole numbers'); end
if any(chr<1 | chr>ct) error('chr must be 1-%d',ct); end

% error checking on pos

if ~isnumeric(pos) error('pos must be numeric'); end
if size(pos,1)==1, pos=pos'; end
if length(pos) ~= length(chr) error('chr and pos arrays must be same length'); end
if any(pos-round(pos)) error('pos must be whole numbers'); end
if any(pos<1) error('pos must be >=1'); end
if any(pos>clen(chr)) error('pos must be within chr length'); end


% computation

offset = [0;cumsum(clen)];
gpos = offset(chr) + pos;
