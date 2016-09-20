function write_hdf5_elements(filename,datasetname,dat,colidx,rowidx)
% WRITE_HDF5_CONT_COLS writes columns of contiguous data over a specified range.
%
%       STAT = GET_HDF5_CONT_COLS(FILENAME,DATASETNAME,DAT,COLRANGE,COLLENGTH)
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
% $Date: 2007-11-29 10:22:55 -0500 (Thu, 29 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$

%% Set the data type to write to file

dtype = class(dat);

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

file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT'); %open with read/write


dataset = H5D.open(file,datasetname);

dataspace = H5D.get_space(dataset);




%%%%%% Check indexing vectors %%%%%%
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


%     %%%%%%%%make the coord array for data hyperslab%%%%%%%%%%%%%%%%%%%%%%%
colcoord = reshape(repmat(colidx',1,length(rowidx))',1,length(rowidx)*length(colidx))-1;
rowcoord = repmat(rowidx,1,length(colidx))-1;

coord = [colcoord; rowcoord];  %make it into a 2XN array


%%%%%%%% select the dataset elements
H5S.select_elements(dataspace,'H5S_SELECT_SET',coord);

% %%%%%%%%%%make memory hyperslab %%%%%%%%%%%%
% rowlength = length(rowidx);
% collength = length(colidx);
% count = [rowlength collength];
%     
% %the mem dataspace
% memspace = H5S.create_simple(2,count,[]);
% 

%% Write file

%H5D.write(dataset,htype,memspace,dataspace,'H5P_DEFAULT',dat);
H5D.write(dataset,htype,'H5S_ALL',dataspace,'H5P_DEFAULT',dat);

%% Close resources

H5S.close(dataspace);
H5S.close(memspace);
H5D.close(dataset);

H5F.close(file)