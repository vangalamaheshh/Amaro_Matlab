function s = remove_trailing_whitespace(s)
s = regexprep(s,'(\s)*$','');
