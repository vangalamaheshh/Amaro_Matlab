function n = editstring_to_n_errors(editstring)
%
% calculates the total number of gaps + mismatches from an editstring
% editstring should be enclosed in brackets
%
% (accepts a cell array of editstrings)
%

if ischar(editstring), editstring = {editstring}; end
nr = length(editstring);
n = nan(nr,1);
dbstop if error
for r=1:nr
  es = editstring{r};
  if strcmp(es,'[-1]'), n(r) = -1; continue; end
  e = sscanf(es(2:end-1),'%d,');
  n(r) = sum(abs([e(2:3:end);e(4:3:end)]));
end
