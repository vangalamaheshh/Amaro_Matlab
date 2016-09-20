function dlist = splitsave_D(D,name,treemap,outpath,segpath,namemap)
%SPLITSAVE_D save D structs for all enumerated subtypes
%
%   
%
%! TODOC

if ~exist('namemap','var') || isempty(namemap)
    namemap = containers.Map;
end

dlist = {}; % list of gcmtypes

if isKey(treemap,name)
    % name is a superclass/interior node => process each child class
    names = treemap(name);
    for i=1:length(names)
        newlist = splitsave_D(D,names{i},treemap,outpath,segpath,namemap);
        dlist = [dlist,newlist];
    end
else
    % name is a TCGA tissue type/leaf node
    dlist = {name};
end
% make list of included gcmtypes into a logical index
keepsamps = false(1,length(D.sis));
for i = 1:length(dlist)
    keepsamps = keepsamps | strcmp({D.sis.gcmtype},dlist{i});
end
% extract samples (if any) and save files
if any(keepsamps)
    % only process leaf nodes and "real" branches (more than one type)
    if ~exist('names','var') || length(names) > 1
        Dsave = reorder_D_cols(D,keepsamps);
        save_D([outpath,'D.',name,'.mat'],Dsave,'-v7.3');
%{
        % write segment files, renaming if necessary
        if isKey(namemap,name)
            segname = namemap(name);
        else
            segname = name;
        end
        write_seg_file([segpath,segname,'.seg'],Dsave);
%}
    end
end
            
