function S = orderfields_last(S,last_flds)
% orders the given fields last; appends the rest of the fields in existing order
%
% Mike Lawrence 2009-11-11

if ischar(last_flds), last_flds = {last_flds}; end

all_flds = fieldnames(S);

if ~isempty(setdiff(last_flds,all_flds)), error('Some of those fields don''t exist'); end

rest_flds = all_flds;
rest_flds(ismember(rest_flds,last_flds)) = [];

S = orderfields(S,[as_column(rest_flds);as_column(last_flds)]);

