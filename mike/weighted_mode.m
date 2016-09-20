function m = weighted_mode(values,weights)
% weighted_mode(values,weights)
%
% returns most frequent value (after weighting)
% in case of a tie, returns smallest of the tied values
%  (matching behavior of native MATLAB "mode" function)
%
% Mike Lawrence 2009-06-25

if length(values)~=length(weights), error('values and weights must be vectors of equal length'); end

[v vi vj] = unique(values);
w = nan(length(v),1);
for i=1:length(v), w(i) = sum(weights(vj==i)); end
mw = max(w);
idx = find(w==mw);
winners = sort(v(idx));
m = winners(1);
