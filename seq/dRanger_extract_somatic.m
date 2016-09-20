function t = dRanger_extract_somatic(t,permissive_flag)
%
% t = dRanger_extract_somatic(t,permissive_flag)
%
% if permissive_flag is set to true,
%   then only polarity-matched evidence from the normal disqualifies.
% if permissive_flag is set to false (default)
%   then ANY evidence from the normal disqualifies.
%
if ~exist('permissive_flag','var'), permissive_flag = false; end

nt = length(t);
nsupp = nan(nt,1);

if permissive_flag
  for i=1:nt, nsupp(i) = length(t{i}.normal_support); end
else
  for i=1:nt
    nsupp(i) = length(t{i}.normal_related_pk{1,1}) +...
               length(t{i}.normal_related_pk{2,1}) +...
               length(t{i}.normal_related_pk{1,2}) +...
               length(t{i}.normal_related_pk{2,2});
  end
end

t = t(nsupp==0);
