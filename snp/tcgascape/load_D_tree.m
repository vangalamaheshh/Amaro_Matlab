function [D markers] = load_D_tree(name,treemap,pathmap,markers,opts)
%LOAD_D_TREE load all the files enumerated from TCGA type hierarchy
%
%   [D MARKERS] = LOAD_D_TREE(NAME,TREEMAP,PATHMAP,MARKERS,OPTS)
%
%! TODOC

D = [];
if isKey(treemap,name)
    % name is a superclass/interior node
    % => load a D for each child class
    names = treemap(name);
    for i=1:length(names)
        [Dnew markers] = load_D_tree(names{i},treemap,pathmap,markers,opts);
        if i==1
            D = Dnew;
        else
            D = merge_Ds(D,Dnew);
        end
    end
else
    % name is a TCGA tissue type/leaf node
    % => load D from segfile (if we can)
    if isKey(pathmap,name)
        segfile = pathmap(name);
        if exist(segfile,'file')
            verbose('Loading CN segments for %s from %s',30,name,segfile);
            D = make_D_from_seg(segfile,markers,opts);
            % use markers in memory for speed
            if ~iscell(markers)
                markers = {D.chrn,D.pos};
            end
            % add TCGA name as sample info gcmtype
            D.sis = struct('gcmtype',cellstr(repmat(name,size(D.dat,2),1)));
        else
            warning('Input segfile %s not found',segfile);
        end
    end
end
            
