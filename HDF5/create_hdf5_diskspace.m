function create_hdf5_diskspace(filename,datasetname,dims,mtype)
% CREATE_HDF5_DISKSPACE creates hdf5 diskspace DATASETNAME in FILENAME for a matrix of MTYPE with size DIMS. 
%
%       CREATE_HDF5_DISKSPACE(FILENAME,DATASETNAME,DIMS,MTYPE)
%
%      ** NOTE: Dims should be passed using the matlab convention ([<col_length> <row_length>]).  To
%      maintain consistency with hdf5write and other files, the created
%      dataspace will be transposed from DIMS.  
%

%       Revisions:
%         19 Dec 07:  Function created (jdobson@broad.mit.edu)
%---
% $Id$
% $Date: 2008-03-21 17:01:30 -0400 (Fri, 21 Mar 2008) $
% $LastChangedBy: jdobson $
% $Rev$


%% Check inputs

if ~ischar(filename) || ~ischar(datasetname)
    error('Inputs FILENAME and DATASETNAME must be char strings');
end

if ~is_full_filename(filename)
    try
    [s,f] = fileattrib(filename);
    filename = f.Name;
    catch
        filename = [pwd '/' filename];
    end
    
end



if size(dims)~= [1 2]
    dims = dims;
    if size(dims)~=[1 2]
        error('DIMS must be a vector of size [1 2]');
    end
end



%% Set the data type


htype = mtype_to_htype(mtype);




%% Define block (select dataspace hyperslab)


    % First, open dataspace
    if exist(filename,'file')
        file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    else
        file = H5F.create(filename,'H5F_ACC_EXCL','H5P_DEFAULT','H5P_DEFAULT'); %fail if file already exists
    end
    
    dims = fliplr(dims);  %convert from matlab style to C-style
    dataspace = H5S.create_simple(2,dims,[]);
    dataset = H5D.create(file,datasetname,htype,dataspace,'H5P_DEFAULT');
    
    H5D.close(dataset);
    H5S.close(dataspace);
    H5F.close(file);
    
    
    
    
    
    

