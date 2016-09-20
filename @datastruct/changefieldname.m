function D = changefieldname(D,oldfield,newfield)
% D = changefieldname(D,oldfield,newfield);

idx = strmatch(oldfield,D.fieldnames,'exact');
D.fieldnames{idx} = newfield;
