function B = expand_slashes_in_genelist(A)
% given list A, takes entries containing "/" (or "//" or "///")
% and expands them to their constituents, returns list B (which may be longer than A)
%
% also removes whitespace
%
% Mike Lawrence 2010-10-07

B = [];
for i=1:length(A)
  a = split(A{i},'/');
  a = regexprep(a,'\s','');
  a(cellfun('isempty',a)) = [];
  B = [B;a];
end
