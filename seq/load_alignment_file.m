function x = load_alignment_file(file,method)

if ~exist('method','var'), method=1; end

switch method

case 1   % 15 seconds

t = load_textfile(file);

idx = find(t=='f');
t(idx) = '+';
t(idx+1) = '1';

idx = find(t=='r');
t(idx) = '-';
t(idx+1) = '1';

d = sscanf(t,'%d');

x = [];
x.chr = d(1:4:end);
x.pos = d(2:4:end);
x.strand = d(3:4:end);
x.id = d(4:4:end);

case 2   % 30 seconds

x = load_struct(file, '%f%f%s%f', 0);
x=rename_field(x,{'col1';'col2';'col3';'col4'},{'chr';'pos';'strand';'id'});

end

x.pos = x.pos+1;
x.start = x.pos;
x.end = x.pos+75;

