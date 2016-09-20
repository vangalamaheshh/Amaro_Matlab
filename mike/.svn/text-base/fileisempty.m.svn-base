function e = fileisempty(fname)
d = dir(fname);
if isempty(d), fprintf('File not found\n'); e=true; return; end
if length(d)>1, fprintf('More than one file matches: using first one\n'); d = d(1); end
e = (d.bytes==0);
