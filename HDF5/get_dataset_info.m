function dinfo = get_dataset_info(fname,dname)
%  GET_DATASET_INFO gets the hdf5 dataset info structure from the hdf5 filename
%  and dataset name.
%
%       DINFO = GET_DATASET_INFO(FNAME,DNAME)
%

%       Revisions:  
%
%           26 Nov 07: Function created by Jen Dobson
%           (jdobson@broad.mit.edu)
%---
% $Id$
% $Date: 2010-04-22 18:36:30 -0400 (Thu, 22 Apr 2010) $
% $LastChangedBy: schum $
% $Rev$


idx = [];

fileinfo = hdf5info(fname);

thislevel = fileinfo.GroupHierarchy;

if isempty(strmatch('/',dname)) % here '/' and not filesep since within HDF5
    dname = ['/' dname];
end

if ~isempty(thislevel.Datasets)
    datasets = thislevel.Datasets;
    idx = strmatch(dname,{datasets.Name},'exact');
end

groups = thislevel.Groups;

while isempty(idx) && ~isempty(groups)

    datasets = groups(1).Datasets;
    idx = strmatch([dname],{datasets.Name});
    if length(groups)>1
        groups = [groups(2:end) groups.Groups];
    else
        groups = groups.Groups;
    end

end
if ~isempty(idx)
    dinfo = datasets(idx);
%! BEGIN patch for one sample schum 2010/4/15
    if dinfo.Rank == 1
        dinfo.Rank = 2;
        dinfo.Dims = [dinfo.Dims 1];
    end
%! END schum hackses instead of doing his taxes ;^)
else
    error('Dataset not found')
end