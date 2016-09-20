function x = cheapperm(n,randseed)
% cheapperm(n,randseed)
%
% quick and dirty version of randperm()
%
% if randseed is provided, initializes twister using that randseed.
%
% Mike Lawrence 2009-05-28

if exist('randseed','var'), rand('twister',randseed); end

s = ceil(sqrt(n));
x = nan(s*s,1);
x(1:n) = 1:n;
x = reshape(x,s,s);
for i=1:7
  p = randperm(s);
  x = x(:,p)';
end
x = reshape(x,1,s*s);
x = x(~isnan(x));
