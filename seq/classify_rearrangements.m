function class = classify_rearrangements(R,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'include_insertions',false);
P = impose_default_value(P,'max_diff_for_insertions',250);
P = impose_default_value(P,'min_diff_for_long_range',1e6);

flds = {'chr1','chr2','pos1','pos2','str1','str2'};
require_fields(R,flds);
R = make_numeric(R,flds);

class = repmat({'unknown'},slength(R),1);

idx = find(R.str1==R.str2);
class(idx) = repmat({'inversion'},length(idx),1);

idx = find(R.str1==0 & R.str2==1);
class(idx) = repmat({'deletion'},length(idx),1);

idx = find(R.str1==1 & R.str2==0);
class(idx) = repmat({'tandem_dup'},length(idx),1);

idx = find(R.pos2-R.pos1>P.min_diff_for_long_range);
class(idx) = repmat({'long_range'},length(idx),1);

if P.include_insertions
  idx = find(R.pos2-R.pos1<P.max_diff_for_insertions);
  class(idx) = repmat({'insertion'},length(idx),1);
end

idx = find(R.chr1~=R.chr2);
class(idx) = repmat({'inter_chr'},length(idx),1);
