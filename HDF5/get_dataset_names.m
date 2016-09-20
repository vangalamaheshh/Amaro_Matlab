function names = get_dataset_names(fname)
%   NAMES = GET_DATASET_NAMES(FNAME) list the datasets in hdf5 file FNAME.

names = [];

fileinfo = hdf5info(fname);

toplevel = fileinfo.GroupHierarchy;

if ~isempty(toplevel.Datasets)
    datasets = toplevel.Datasets;

    names = {datasets.Name};
end

groups = toplevel.Groups;

while ~isempty(groups)
    names = [names {groups(1).Datasets.Name}];
    if length(groups)>1
        groups = [groups(2:end) groups(1).Groups];
    else
        groups = [groups(1).Groups];
    end
end
