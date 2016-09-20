function Q = parse_qltout(file)
X = load_textfile(file);
cr = find(X==char(10));
nl = length(cr);

A = cell(nl,11);

matchstring = '^QUERY\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(\d*)\t(.*)$';

pos = 1;
n=1;
for l=1:nl
  L = X(pos:cr(l)-1); pos = cr(l)+1;
  if strncmp(L,'QUERY',5)
    tmp = regexp(L,matchstring,'tokens');
    if isempty(tmp{1}), fprintf('Error parsing string\n'); keyboard; end
    A(n,1:10) = tmp{1};
    A{n,10} = ['[' regexprep(A{n,10},'\t',',') ']'];
    A{n,11} = '';
    n=n+1;
  elseif strncmp(L,'#',1)
    A{n-1,11} = L;
  else
    fprintf('Error parsing line\n'); keyboard;
  end
end

Q=[];
Q.id = str2double(A(:,1));
Q.qstart = str2double(A(:,2)) + 1;
Q.qend = str2double(A(:,3));
Q.length = str2double(A(:,4));
Q.rc = str2double(A(:,5));
Q.chr = str2double(A(:,6));
Q.start = str2double(A(:,7)) + 1;
Q.end = str2double(A(:,8));
Q.chrlength = str2double(A(:,9));
Q.editstring = A(:,10);
Q.comment = A(:,11);

Q = reorder_struct(Q,1:n-1);
