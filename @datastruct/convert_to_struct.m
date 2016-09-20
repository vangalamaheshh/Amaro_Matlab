function S = convert_to_struct(D)



fields = fieldnames(D);

memtype = D.storagetype;
isdisk = strmatch('disk',memtype);
D = convert_to_memfield(D,fields(isdisk));
    

for k = 1:length(fields)
    S.(fields{k}) = D.fielddata{k};
end

