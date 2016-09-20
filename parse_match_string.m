function x=parse_match_string(match)
x=strsplit(match,'>');
x=strsplit(x{2},'views<');
x=str2double(regexprep(x{1},',',''));
end