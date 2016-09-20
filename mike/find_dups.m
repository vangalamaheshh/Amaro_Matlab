function [dup idx ct] = find_dups(varargin)

[u tmp uj] = unique_combos(varargin{:});
h = histc(uj,1:length(u));
z = find(h>1);
dup = u(z);
ct = h(z);
idx = cell(length(z),1);
for i=1:length(z)
  idx{i} = find(uj==z(i));
end
