function write_hdf5_block(filename,datasetname,dat,colidx,rowidx)
% GET_HDF5_BLOCK returns data indexed by indices COLIDX and ROWIDX when
% indices define a block of data in hdf5 file.
%
%       WRITE_HDF5_BLOCK(FILENAME,DATASETNAME,DAT,COLIDX,ROWIDX)
%
%       GET_HDF5_BLOCK returns the elements of the dataset in file FILENAME
%       and dataset DATASETNAME defined by COLIDX and ROWIDX.  DTYPE is the
%       data type of the dataset.  GET_HDF5_BLOCK should be used when
%       accessing a contiguous block of data or accessing partial columns
%       or rows of data.  If a contiguous block or partial rows/columns
%       cannot be found, GET_HDF5_ELEMENTS is called instead.
%
%

%       Revisions:
%           28 Nov 07: Function created by Jen Dobson
%           (jdobson@broad.mit.edu).
%
%           4 Dec 07:  Changed memspace selection to specify size using a
%           block (last argument) rather than count (second-to-last
%           argument) in:
%               H5S.select_hyperslab(memspace,'H5S_SELECT_SET',[0 0],[],[1 1],count);
%           11 Apr 11: substantial rewrite using memory buffers for indexing 
%           instead of using HS5_SELECT_OR to create an HDF5 index (schum)
%---
% $Id$
% $Date: 2011-04-11 15:18:29 -0400 (Mon, 11 Apr 2011) $
% $LastChangedBy: schum $
% $Rev$


%% check indexing vectors

if size(colidx,1) > size(colidx,2)
    colidx = colidx';
end
if size(rowidx,1) > size(rowidx,2)
    rowidx = rowidx';
end
if size(rowidx,1) ~= 1
    error('ROWIDX must be a vector');
end
if size(colidx,1) ~= 1
    error('COLIDX must be a vector');
end


%% set the data type to write to file

dtype = class(dat);
htype = mtype_to_htype(dtype);

%% open dataset

file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT'); %open with read/write access
dataset = H5D.open(file,datasetname);
dataspace = H5D.get_space(dataset);

%% test inputs for conditions for optimizations

% find index ranges (block limits)
maxr = max(rowidx);
minr = min(rowidx);
maxc = max(colidx);
minc = min(colidx);

% if entire block will be filled, no need to read before writing
rows_hit = false(1,maxr-minr+1);
rows_hit(rowidx-minr+1) = true;
cols_hit = false(1,maxc-minc+1);
cols_hit(colidx-minc+1) = true;
mapping_is_onto = all(rows_hit) & all(cols_hit);

% if indexing is identity mapping, no need to reorder data
contiguous_indices = isequal(rowidx,rowidx(1):rowidx(end)) & ...
                     isequal(colidx,colidx(1):colidx(end));


%% read current disk image into memory
%! TODO (optimize) only if not contiguous

% select a hyperlab that spans the range of the indices
offset = [minc-1,minr-1]; 
block = [maxc maxr] - offset;
count = [1 1];
stride = [];
H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,stride,count,block);

% create the memory space
memspace = H5S.create_simple(2,block,[]);
% transfer data from disk to memory buffer if we need to

%% write the data
if contiguous_indices
    % contiguous block: simple write
    H5D.write(dataset,htype,memspace,dataspace,'H5P_DEFAULT',dat);
else
    % if there are gaps in the indices, need to read data to fill in
    if ~mapping_is_onto
        diskdat = H5D.read(dataset,htype,memspace,dataspace,'H5P_DEFAULT');
    end
    % reorder data and write
    diskdat(rowidx-offset(2),colidx-offset(1)) = dat;
    H5D.write(dataset,htype,memspace,dataspace,'H5P_DEFAULT',diskdat);
end

%% Close resources
H5S.close(dataspace);
H5D.close(dataset);
H5F.close(file)
