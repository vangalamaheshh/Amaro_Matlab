function S = fillblanks(S,fill)
%
% fillblanks(S,fill)
%
% S is a cell array of strings
% fillblanks replaces all empty strings with the text given as "fill"
%

if ~isempty(S)
  idx = find(cellfun('isempty',S));
  S(idx) = repmat({fill},length(idx),1);
end
