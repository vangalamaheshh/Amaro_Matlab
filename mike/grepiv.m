function keep = grepiv(pattern,strings,flag)
if ~exist('flag','var'), flag=0; end
keep = grepvi(pattern,strings,flag);
