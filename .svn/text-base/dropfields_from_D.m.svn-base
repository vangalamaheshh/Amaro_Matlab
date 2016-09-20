function D = dropfields_from_D(D,fields_to_drop)

if ~exist('fields_to_drop','var')
    fields_to_drop = {'history','orig','origidx','gorigidx'};
end
fields_to_drop = fields_to_drop(isfield(D,fields_to_drop));
if ~isempty(fields_to_drop)
    D = rmfield(D,fields_to_drop);
end