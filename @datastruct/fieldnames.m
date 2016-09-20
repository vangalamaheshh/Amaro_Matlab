function fn = fieldnames(ds)
%FIELDNAMES list the fieldnames of datastruct DS.
%
%   FN = FIELDNAMES(DS)

fn = ds.fieldnames';

if size(fn,2) > size(fn,1)
    fn = fn';
end

