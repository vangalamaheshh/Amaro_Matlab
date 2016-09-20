function X = load_boulder_file(fname)
% Mike Lawrence 2010-04-06

x = load_lines(fname);
div = strcmp(x,'=');
rec = 1+cumsum(div);
x(div) = []; rec(div) = [];
t = parse(x,'(.*)=(.*)',{'key';'val'});

[u ui uj] = unique(t.key);
q = cell(rec(end),length(u));
for i=1:slength(t), q{rec(i),uj(i)} = t.val{i}; end

X = [];
flds = genfieldname(u);
[tmp ord] = sort(ui);
for j=1:length(u), X = setfield(X,flds{ord(j)},q(:,ord(j))); end
