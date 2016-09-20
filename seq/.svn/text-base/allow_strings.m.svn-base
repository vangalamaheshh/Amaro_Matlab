function list = allow_strings(list, allowed, varargin)
%
% allow_strings(list, allowed, parameters)
%
% Ensures only the strings in allowed appear in list.
%
% If 'case_insensitive' appears among parameters,
%    then case is ignored in making the comparisons,
%    but the list is updated to have the case that appears in allowed.
%
% Mike Lawrence 2008-05-27
%

case_insensitive = false;

for i=1:length(varargin)
  if strcmpi(varargin{i}, 'case_insensitive')
    case_insensitive = true;
  else
    error('Unrecognized option %s', varargin{i});
  end
end

for i=1:length(list)
  string_ok = false;
  for j=1:length(allowed)
    if case_insensitive
      match = strcmpi(allowed{j}, list{i});
    else
      match = strcmp(allowed{j}, list{i});
    end
    if match
      list{i} = allowed{j};
      string_ok = true;
      break;
    end
  end
  if ~string_ok
    keyboard
    error('Non-allowed string %s', list{i});
  end
end
