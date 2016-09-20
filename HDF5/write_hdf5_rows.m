function stat = write_hdf5_cont_rows(filename,datasetname,dat,rowidx,rowlength,dtype)
% WRITE_HDF5_CONT_ROWS writes rows of data indexed by input argument ROWIDX.
%
%       STAT = WRITE_HDF5_CONT_ROWS(FILENAME,DATASETNAME,DAT,ROWIDX,ROWLENGTH,DTYPE)
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
% $Date: 2007-11-20 14:40:15 -0500 (Tue, 20 Nov 2007) $
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


count = [rowlength length(rowidx)];

onerowct = [rowlength 1];

%initialize the dataset's hyperslab
offset = [0 rowidx(1)-1];
H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,[],onerowct,[]);
%append to the dataset's hyperslab
for k = rowidx(2:end)
    offset = [0 k-1];
    H5S.select_hyperslab(dataspace,'H5S_SELECT_OR',offset,[],onerowct,[]);
end
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
        