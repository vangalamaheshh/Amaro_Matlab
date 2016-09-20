function dlist = write_igv_segs(D,name,treemap,segpath,namemap,minsamples)
% WRITE_IGV_SEGS write segfiles for all enumerated subtypes
%
%   DLIST = write_igv_segs(D,NAME,TREEMAP,SEGPATH,NAMEMAP,MINSAMPLES)
%
%

if ~exist('namemap','var') || isempty(namemap)
    namemap = containers.Map;
end

if ~exist('minsamples','var') || isempty(minsamples)
    minsamples = 40;
end

dlist = {}; % list of gcmtypes

if isKey(treemap,name)
    % name is a superclass/interior node => process each child class
    names = treemap(name);
    for i=1:length(names)
        newlist = write_igv_segs(D,names{i},treemap,segpath,namemap);
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
if sum(keepsamps) >= minsamples
    % only process leaf nodes and "real" branches (more than one type)
    if ~exist('names','var') || length(names) > 1
        Dsave = reorder_D_cols(D,keepsamps);        
        % write segment files, renaming if necessary
        if isKey(namemap,name)
            segname = namemap(name);
        else
            segname = name;
        end
        fname = regexprep(lower(segname),' ','_');
        write_seg_file([segpath,fname,'.seg'],Dsave);
    end
end
            
