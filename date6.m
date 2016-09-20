function str6 = date6(timestamp)
%DATE6 create 6-digit date string for current time
if ~exist('timestamp','var') || isempty(timestamp)
    timestamp = now;
end
str8 = date8(timestamp);
str6 = str8(3:end);
