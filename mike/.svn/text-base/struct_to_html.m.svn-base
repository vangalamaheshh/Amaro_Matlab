function H = struct_to_html(L,P)
% H = struct_to_html(L,P)
%
% Mike Lawrence 2010-04-26

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'flds1',fieldnames(L));
P=impose_default_value(P,'flds2',P.flds1);
P=impose_default_value(P,'maxrows',slength(L));
P=impose_default_value(P,'table_class','its');

if ~isnumeric(P.maxrows), P.maxrows = str2double(P.maxrows); end
nf = length(P.flds1);
nl = min(slength(L),P.maxrows);

P=impose_default_value(P,'widths',50*ones(1,nf));
P=impose_default_value(P,'formnum',zeros(1,nf));
P=impose_default_value(P,'sigfigs',3*ones(1,nf));  % used if formnum==1
P=impose_default_value(P,'numwidth',8*ones(1,nf)); % used if formnum==1
P=impose_default_value(P,'lessthan',false(nl,nf)); % used if formnum==1
P=impose_default_value(P,'decplc',zeros(1,nf));    % used if formnum==0
P=impose_default_value(P,'commas',zeros(1,nf));    % used if formnum==0
P=impose_default_value(P,'justfy',-ones(1,nf));
P=impose_default_value(P,'rowcolor',cell(nl,1));
P=impose_default_value(P,'columntextcolor',cell(nf,1));

if isempty(P.table_class)
  H = '<table>';
else
  H = ['<table class="' P.table_class '">'];
end

H = [H '<tr>'];

X = cell(nl,nf);
for f=1:nf
  ff = ['%0.' num2str(P.decplc(f)) 'f'];
  tmp = getfield(L,P.flds1{f});
  if isnumeric(tmp)
    for i=1:nl
      if isnan(tmp(i))
        X{i,f} = '&mdash;';
      else
        if P.formnum(f)
          X{i,f} = format_number(tmp(i),P.sigfigs(f),P.numwidth(f));
          if P.lessthan(i,f)
            X{i,f} = ['&lt;' X{i,f}];
          end
        else
          X{i,f} = sprintf(ff,tmp(i));
          if P.commas(f)==1, X{i,f} = comma_format(X{i,f}); end
        end
      end
    end
  else
    X(:,f) = tmp(1:nl);
  end
end

for f=1:nf
  H = [H '<th'];
  if P.justfy(f)==1, H = [H ' align=right']; end
  if P.justfy(f)==-1, H = [H ' align=left']; end
  H = [H ' width="' num2str(P.widths(f)) '">'];
  if ~isempty(P.columntextcolor{f}), H = [H '<font color=' P.columntextcolor{f} '>']; end
  H = [H P.flds2{f}];
  if ~isempty(P.columntextcolor{f}), H = [H '</font>']; end
  H = [H '</th>'];
end
H = [H '</tr>'];

%fidx = find(strcmp(P.flds1,'fusion'));
Rows = cell(nl,1);
for i=1:nl, if ~mod(i,1000), fprintf('%d/%d ',i,nl); end
  R = '<tr>';
  for f=1:nf
    R = [R '<td'];
    if ~isempty(P.rowcolor{i}), R = [R ' bgcolor=' P.rowcolor{i}]; end
    if P.justfy(f)==1, R = [R ' align=right']; end
    if P.justfy(f)==-1, R = [R ' align=left']; end
    R = [R '>'];
%    redflag = (f==fidx && L.fusion_priority(i)>=1);
%    if redflag, R = [R '<font color="red"><b>']; end
    if ~isempty(P.columntextcolor{f}), R = [R '<font color=' P.columntextcolor{f} '>']; end
    R = [R X{i,f}];
%    if redflag, R = [R '</b></font>']; end
    if ~isempty(P.columntextcolor{f}), R = [R '</font>']; end
    R = [R '</td>'];
  end
  R = [R '</tr>'];
  Rows{i} = R;
end, if i>=1000, fprintf('\n'); end

H = [H cat(2,Rows{:}) '</table>'];
