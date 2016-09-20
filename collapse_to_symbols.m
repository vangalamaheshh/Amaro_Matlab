function X=collapse_to_symbols(X,collapse_type)

if ~exist('collapse_type','var')
  collapse_type='median';
end

% remove rows w/o symbols
X=reorder_D_rows(X,find(cellfun('isempty',regexp(X.gsymb,'---'))));

% collapse by symbol
[us,ui,uj]=unique(strvcat(lower(X.gsymb)),'rows');
[suj,suji]=sort(uj);
X=reorder_D_rows(X,suji);
X=add_D_sup(X,'SYM','Symbol id',suj','rows');
X=collapse_D(X,1,collapse_type);

X.gacc=X.gsymb;

