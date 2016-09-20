function [a]=chrom2num(chr)
 c=chr;
 if (length(c)<1), a=[]; return; end
 cx = regexprep(cellstr(c), 'X', '23', 'ignorecase');
 cx = regexprep(cx, 'Y', '24', 'ignorecase');
 cx = regexprep(cx, 'MT', '25', 'ignorecase');
 cx = regexprep(cx, 'M', '25', 'ignorecase');
 a=str2num(char(cx))';
 a=a(:);
end