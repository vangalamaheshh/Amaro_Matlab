function s = combined_std(mu,sd,n)
% combined_std(mu,se,n)
%
% given the means, standard deviations, and n's of a set of groups,
% computes the standard deviation of the combined group
% (divides by N-1)
%
% Mike Lawrence 2008-08-13

if length(mu)~=length(sd) || length(mu)~=length(n)
  error('input vectors must be of equal length');
end

sxm2 = (sd.^2) .* (n-1);   % = sum(x-mu)^2
sx2 = sxm2 + (n .* mu.^2); % sum(x^2)
ssx2 = sum(sx2);
sn = sum(n);
totmu = sum(n.*mu)/sn;

s2 = (ssx2 - sn*totmu^2) / (sn-1);
s = sqrt(s2);

