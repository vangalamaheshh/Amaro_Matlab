function x = prefix(x,pref)

if ~iscellstr(x), error('expects cellstr for first argument'); end
if ~ischar(pref), error('expects char for second argument'); end

x = regexprep(x,'^(.*)$',[pref '$1']);
