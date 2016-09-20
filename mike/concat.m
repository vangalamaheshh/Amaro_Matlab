function c = concat(strings,separator)
%
% concat(strings,separator)
%
% joins strings into one string, separated with the specified character/string
%
% Mike Lawrence 2008-05-01

c='';
for i=1:length(strings)
  if ischar(strings(i))
    c=[c strings(i)];
  elseif isnumeric(strings(i))
    c=[c num2str(strings(i))];
  elseif isnumeric(strings{i})
    c=[c num2str(strings{i})];
  elseif iscell(strings(i))
    c=[c strings{i}];
  else
    error('concat does not know how to handle that kind of input');
  end
  if i<length(strings)
    c=[c separator];
  end
end
