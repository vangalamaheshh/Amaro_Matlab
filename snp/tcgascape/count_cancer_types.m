function type2count = count_cancer_types(name,sample_types,treegen,type2count)
%
%! TODOC

if ~exist('type2count','var') || isempty(type2count)
    type2count = containers.Map;
end

if isKey(treegen,name)
    % name is a superclass/interior node => process each child class
    subtypes = treegen(name);
    count = 0;
    for i=1:length(subtypes)
        type2count = count_cancer_types(subtypes{i},sample_types,treegen,type2count);
        count = count + type2count(subtypes{i});
    end
    type2count(name) = count;
else
    % name is a TCGA tissue type/leaf node
    type2count(name) = sum(strcmp(sample_types,name));
end

            
