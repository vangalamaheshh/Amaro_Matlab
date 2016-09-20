function D = addfield(D,field,dims,type,dm,hdf5dir)


%get hdf5 directory of first diskfield
if ~exist('type','var') || isempty(type)
    type = 'double';
end

if strcmp(D.fieldnames,field)
    error('Field already exists');
end

if ~exist('dm','var') || isempty(dm)  %if storage type is not specified, set it to memory unless dims matches another diskfield's dims
    dm = 'memory';
    for k = strmatch('disk',D.storagetype)'
        if dims == D.fielddata{k}.Dims
            dm = 'disk';
            continue;
        end
    end
end

if strcmpi('disk',dm)

    if ~exist('hdf5dir','var') || isempty(hdf5dir)
        idx = strmatch('disk',D.storagetype);
        f = D.fieldnames{idx(1)};
        filename = regexprep(get_datafile(D,f),f,field);
    else
        if ~strcmp(hdf5dir(end),filesep)
            hdf5dir = [hdf5dir filesep];
        end
        filename = [hdf5dir field '.h5'];
        if ~is_full_filename(filename)
            filename = [pwd filesep filename];
        end
    end
    D = add_diskfield(D,filename,field,dims,type);
    
else
    data = cast(zeros(dims),type);
    D = add_memfield(D,field,data);
end

