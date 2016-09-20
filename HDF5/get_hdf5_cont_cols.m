function dat = get_hdf5_cont_cols(filename,datasetname,colrange,collength,dtype)
% GET_HDF5_CONT_COLS returns the columns of contiguous data over a specified range.
%
%       DAT = GET_HDF5_CONT_COLS(FILENAME,DATASETNAME,COLRANGE,COLLENGTH,DTYPE)
%
%       Note: When referencing data using the hdf5 C library, use the C convention of
%       specifying row dimension first, then column dimension.  
%
%      

%       Revisions:
%           20 Nov 07: Function created by Jen Dobson
%           (jdobson@broad.mit.edu).
%
%           27 Nov 07: Cleaned up dataspace and memspace definitions:
%           Changed definition of dataspace hyperslab (specified block of
%           data rather than a count);  Write to whole memory space rather
%           than defining a hyperslab in memory.
%
%---
% $Id$
% $Date: 2007-11-29 10:22:55 -0500 (Thu, 29 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$


%% Set the data type to write to file

switch dtype
    case 'single'
        htype = 'H5T_NATIVE_FLOAT';
    case 'double'
        htype = 'H5T_NATIVE_DOUBLE';
    case 'char'
        htype = 'H5T_NATIVE_CHAR';
    case 'logical'
        htype = 'H5T_NATIVE_HBOOL';
    case 'int16'
        htype = 'H5T_NATIVE_SHORT';
    case 'uint16'
        htype = 'H5T_NATIVE_USHORT';
    case 'int32'
        htype = 'H5T_NATIVE_INT';
    case 'uint32'
        htype = 'H5T_NATIVE_UINT';
    case 'int64'
        htype = 'H5T_NATIVE_LONG';
    case 'uint64'
        htype = 'H5T_NATIVE_ULONG';
    otherwise
        htype = dtype; %not a matlab class
end


%% Open file and create data spaces

file = H5F.open(filename,'H5F_ACC_RDONLY','H5P_DEFAULT'); %open with read only access


dataset = H5D.open(file,datasetname);

dataspace = H5D.get_space(dataset);


if isscalar(colrange)
    offset = [colrange-1 0];
    colblock = [1 collength];
elseif size(colrange) == [1 2]
    startidx = colrange(1);
    endidx = colrange(2);
    if startidx > endidx
        error('Second element of 1x2 input vector colrange must be greater than first element')
    end
    offset = [startidx-1 0];
    colblock = [endidx-startidx+1 collength];
else
    error('Input parameter COLRANGE must be scalar or 1x2 row vector')
end

    


%the dataset's hyperslab
H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,[],[1 1],colblock);
    %stride and block are left as []

    
%%%%%%%%%%make memory hyperslab %%%%%%%%%%%%

%the mem dataspace
memspace = H5S.create_simple(2,colblock,[]);

%% Read data

dat = H5D.read(dataset,htype,memspace,dataspace,'H5P_DEFAULT');

%% Close resources


H5S.close(dataspace);
H5S.close(memspace);
H5D.close(dataset);

H5F.close(file)
  