function stat = write_hdf5_cont_rows(filename,datasetname,dat,rowrange,rowlength,dtype)
% WRITE_HDF5_CONT_ROWS writes rows of contiguous data over a specified range.
%
%       STAT = WRITE_HDF5_CONT_ROWS(FILENAME,DATASETNAME,DAT,ROWRANGE,ROWLENGTH,DTYPE)
%
%       Note: When referencing data using the hdf5 C library, use the C convention of
%       specifying row dimension first, then column dimension.  
%
%

%       Revisions:
%           20 Nov 07: Function created by Jen Dobson
%           (jdobson@broad.mit.edu).
%---
% $Id$
% $Date: 2007-11-20 12:47:33 -0500 (Tue, 20 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$


%% If dtype is given, check class of data
if exist('dtype','var')
    if ~strcmp(class(dat),dtype)
        error('Class of input variable DAT does not match input variable DTYPE')
    end
else
    dtype = class(dat);
end

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
        error('Input variable DAT: Unrecognized class.')
end

%% Open file and create data spaces

file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT'); %open with read only access


dataset = H5D.open(file,datasetname);

dataspace = H5D.get_space(dataset);



if isscalar(rowrange)
    offset = [0 rowrange-1];
    count = [rowlength 1];
elseif size(rowrange) == [1 2]
    startidx = rowrange(1);
    endidx = rowrange(2);
    if startidx > endidx
        error('Second element of 1x2 input vector rowrange must be greater than first element')
    end
    offset = [0 startidx-1];
    count = [rowlength endidx-startidx+1];
else
    error('Input parameter ROWRANGE must be scalar or 1x2 row vector')
end


%the dataset's hyperslab
H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,[],count,[]);
    %stride and block are left as []
    
    
%the mem dataspace
memspace = H5S.create_simple(2,count,[]);


%the mem hyperslab
H5S.select_hyperslab(memspace,'H5S_SELECT_SET',[0 0],[],count,[]);


%% Write file
H5D.write(dataset,htype,memspace,dataspace,'H5P_DEFAULT',dat);

%% Close resources

H5S.close(dataspace);
H5S.close(memspace);
H5D.close(dataset);

H5F.close(file)
        