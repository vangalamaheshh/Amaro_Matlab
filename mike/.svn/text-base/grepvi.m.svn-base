function keep = grepvi(pattern,strings,flag)
% grepvi(pattern,strings,flag)
%
% case-insensitive inverse grep
%
% Mike Lawrence 2010-02-04

if ~exist('flag','var'), flag=0; end

toss = grep(upper(pattern),upper(strings),1);
keep = setdiff((1:length(strings))',toss);
if ~flag, keep = strings(keep); end
