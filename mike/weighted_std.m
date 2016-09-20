function s = weighted_std(values,weights)
% weighted_std(values,weights)
%
% standard deviation of weighted values
% (divides by N-1, where weights are assumed to be
% equal to number of occurrences)
%
% Mike Lawrence 2008-06-21

if length(values)~=length(weights)
  error('values and weights must be vectors of equal length');
end

if length(values)==1, s=0; return; end

m = weighted_mean(values,weights);
d2 = (values-m).^2;
t = d2 .* weights;
v = sum(t) / (sum(weights)-1);
s = sqrt(v);

