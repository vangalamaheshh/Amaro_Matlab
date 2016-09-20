function D = addfield(D,field,dims,type,dm,hdf5dir)

if exist('dm','var') && ~isempty('dm')
    warning('Storage type is automatically memory for data in struct objects')
end

if exist('hdf5dir','var') && ~isempty('hdf5dir')
    warning('HDF5 directory input argument ignored for struct object')
end


D.(field) = cast(zeros(dims),type);