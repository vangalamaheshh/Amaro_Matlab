function addtree(rootdir)
% ADDTREE add directory and all subdirectories to the path.
if ~isdeployed
    p = genpath(rootdir);
    p = textscan(p,'%[^:]','Delimiter',':');
    bool = cellfun(@isempty,(regexp(p{:},'.svn')));
    p = p{1}(bool);
    p = strcat(p,':');
    addpath(cell2mat(p'));
end