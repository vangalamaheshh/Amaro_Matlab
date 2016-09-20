function s = ismemfield(D,fields)
%ISDISKFIELD returns true for disk fields.
%
%   S = ISDISKFIELD(D,FIELDS) for datastruct D and cell array of fields of
%   D FIELDS, returns S, an array of ones and zeros with a one at locations where FIELDS is a disk field and zero where FIELDS is a mem field.




if ~iscell(fields)
    fields = {fields};
end


if ~strcmp(class(D),'datastruct')
    error('D is not of class DATASTRUCT')
end



if length(intersect(fieldnames(D),fields)) ~= length(unique(fields))
    error('One or more fields is not a field of input datastruct')
end


[int,dum,idx] = intersect(fields,D.fieldnames);

s = strcmp(D.storagetype{idx},'memory');


