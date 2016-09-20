function isleafmap = leaf_cancer_types(name,treegen,isleafmap)
%LEAF_CANCER_TYPES enumerate the "independent" (leaf) cancer types
%
%   ISLEAFMAP = leaf_cancer_types(NAME,TREEGEN,ISLEAFMAP)
%
% Returns a containers.Map with a key for each leaf subtypes under
% the NAME, for the type hierarchy defined by the containers.Map TREEGEN.
% If optional ISLEAFMAP container.Map is provided, leaf types are added to
% ISLEAFMAP instead of being created.

if ~exist('isleafmap','var') || isempty(isleafmap)
    isleafmap = containers.Map;
end

if isKey(treegen,name)
    % name is a superclass/interior node
    isleafmap(name) = false;
    % process each child class
    subtypes = treegen(name);
    for i=1:length(subtypes)
        isleafmap = leaf_cancer_types(subtypes{i},treegen,isleafmap);
    end
else
    % name is a TCGA tissue type/leaf node
    isleafmap(name) = true;
end

            
