function m = understand_nested_target_sets(T,cutoff)
% understand_nested_target_sets(T,cutoff)
%
% Input should be a cell array of target-sets,
%   each member of which is a struct
%   that includes the following fields: chr, start, end;
%   Target-sets should be ordered by their nesting,
%   with the largest (outermost) set first.
%
% Output is an array of integers,
%   one for each target in the largest target set,
%   telling the smallest (innermost) set each target belongs to.
%
% "Belongs to" = at least "cutoff" fraction (default=0.90)
%   of the bases in the target are also covered jointly by
%   the targets of the smaller set.
%
% For cases where a target belongs to set1 and set3 but not set2
%   (i.e. the target-sets are not TRULY nested),
%   the behavior of the algorithm is to return "3" nonetheless.
%
% Mike Lawrence 2009-06-11

if ~exist('cutoff','var'), cutoff = 0.9; end

if ~iscell(T), error('input should be a cell array'); end
ns = length(T);
if ns<2, error('At least two target sets are required'); end
for s=1:ns
  if ~isstruct(T{s}), error('input should be a cell array of structs'); end
  require_fields(T{s},{'chr','start','end'});
  if isnan(slength(T{s})), error('structs should have uniform length'); end
  if ~isnumeric(T{s}.chr), error('chr should be numeric'); end
  if ~isnumeric(T{s}.start), error('start should be numeric'); end
  if ~isnumeric(T{s}.end), error('end should be numeric'); end
  if any(T{s}.chr<1 | T{s}.chr>24 | isnan(T{s}.chr)) error('chr should be 1-24'); end
end

chrlen = load_chrlen;

nt = slength(T{1});   % number of targets in outermost set
m = ones(nt,1);
for s=2:ns
  fprintf('EXAMINING TARGET SET %d/%d\n',s,ns);
  for chr=1:24
    fprintf('  chr%d\n',chr);
    cov = false(chrlen(chr),1);
    cidx = find(T{s}.chr==chr);
    for j=1:length(cidx); i=cidx(j);
      st = max(1,T{s}.start(i));
      en = min(chrlen(chr),T{s}.end(i));
      cov(st:en) = true;
    end
    cidx = find(T{1}.chr==chr);
    for j=1:length(cidx); i=cidx(j);
      st = max(1,T{1}.start(i));
      en = min(chrlen(chr),T{1}.end(i));
      fcov = mean(double(cov(st:en)));
      if fcov>=cutoff, m(i) = s; end
    end
  end
end
