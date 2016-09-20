function S = load_cnseg_file(fname)
% Mike Lawrence 2009-10-16

if ~exist(fname,'file'), error('Not found: %s',fname); end

f = fopen(fname);
l = fgetl(f);
fclose(f);
tabs = find(l==char(9));
numcols = length(tabs)+1;
if numcols~=6, error('Invalid seg file'); end

l2 = l(tabs(1)+1:end);
d = sscanf(l2,'%d');
if length(d)==5, header_lines = 0; else header_lines = 1; end

x = read_table(fname,'%s%f%f%f%f%f',char(9),header_lines);

S = [];
S.id = x.dat{1};
S.chr = x.dat{2};
S.start = x.dat{3};
S.end = x.dat{4};
S.nprobes = x.dat{5};
S.segmean = x.dat{6};

S.ratio = 2.^S.segmean;
