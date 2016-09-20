function s = bp2str(b,decimal_places)
%
% given a count of basepairs,
% returns a string like '2bp', '2Kb', '2Mb'
%
% Mike Lawrence 2009-02-10

if ~exist('decimal_places','var'), decimal_places = 0; end

fs = ['%0.' num2str(decimal_places) 'f'];

x{1} = sprintf([fs 'bp'],b);
if b>=1000, x{2} = sprintf([fs 'Kb'],b/1000); end
if b>=1000000, x{3} = sprintf([fs 'Mb'],b/1000000); end

l = cellfun('length',x);
[tmp i] = min(l);
s = x{i};

