function D=read_mit_sin_file(D,fname,cols)

% Cel_file_name NAME col1 col2 ....
try
  if strcmp(lower(fname((end-2):end)),'xls')
    [n,t,r]=xlsread(fname);
  else
    r=read_dlm_file(fname);
    r=cat(1,r{:});
    r=regexprep(r,'"','');
  end
catch
end
disp('Case Insensitive!');
r=lower(r);

[Mt,m1,m2]=match_string_sets(lower(D.sdesc),r(2:end,strmatch('name',r(1,:),'exact')));

if length(m1)~=size(D.dat,2)
  disp([ 'No match for:']);
  catn(D.sdesc(setdiff(1:size(D.dat,2),m1),:));
end
if length(m2)~=(size(r,1)-1)
  disp([ 'Unmatched in SIN file']);
  catn(r(1+setdiff(1:(size(r,1)-1),m2),strmatch('name',r(1,:),'exact')));
end

D.scan=cell(size(D.dat,2),1);
D.scan(m1)=r(1+m2,strmatch('cel_file_name',r(1,:),'exact'));

if ischar(cols)
  cols={cols};
end

if iscell(cols)
  [Nt,n1,n2]=match_string_sets(r(1,3:end),lower(cols));
  assert(length(n2)==length(cols));
  disp([ 'Matched: ']);
  disp([ cellstr(num2str(as_column(n1))) cols ]);
  cols=n1;
end

for i=as_row(cols)
  type=r{1,i+2};
  [c,ci,cj]=unique_keepord(strvcat(r{2:end,i+2}),'rows');
  tmp=[ mat2cell((1:size(ci,1))',ones(size(ci,1),1),1) r(ci+1,i+2)]';
  desc=[ type ': ' sprintf('%d-%s/',tmp{:})];
  desc=desc(1:(end-1)); % remove last /
  v=nan(1,size(D.dat,2));
  v(m1)=cj(m2);
  D=add_D_sup(D,desc,desc,v,'col');
end

