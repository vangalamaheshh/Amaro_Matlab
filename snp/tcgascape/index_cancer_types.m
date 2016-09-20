function type2index = index_cancer_types(name,sample_types,treegen,type2index)
%SPLITSAVE_D save D structs for all enumerated subtypes
%
%! TODOC

if ~exist('type2index','var') || isempty(type2index)
    type2index = containers.Map;
end

if isKey(treegen,name)
    % name is a superclass/interior node => process each child class
    subtypes = treegen(name);
    index = false(size(sample_types));
    for i=1:length(subtypes)
        type2index = index_cancer_types(subtypes{i},sample_types,treegen,type2index);
        index = index | type2index(subtypes{i});
    end
    type2index(name) = index;
else
    % name is a TCGA tissue type/leaf node
    type2index(name) = strcmp(sample_types,name);
end

            
