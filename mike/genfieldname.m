function b = genfieldname(a)
% genfieldname(a)
%
% similar to native MATLAB function genvarname,
%   but leaves keywords alone.
%
% Mike Lawrence 2009-02-04

b = genvarname(a);

if iscell(b)
  for i=1:length(b), if iskeyword(a(i)), b(i)=a(i); end, end
else
  if iskeyword(b), b=a; end
end
