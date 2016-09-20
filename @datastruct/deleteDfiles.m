function deleteDfiles(d,field)
% SUC = deleteDfiles(D,FIELD)
%  Delete HDF5 files associated with datastruct.
%
%   If FIELD argument is passed, only files associated with that field are
%   deleted.
%
%  Recommended use is after making a temporary copy of a datastruct with
%  copyD.  Use deleteDfiles when done with the copy before clearing the
%  copy.



if iswriteprotected(d)
    error('Cannot delete hdf5 files associated with datastruct object.  Datastruct object is write protected.')
end


if ~exist('field','var') || isempty(field)

    didx = strmatch('disk',d.storagetype);


    diskfielddata = {d.fielddata{didx}};

    % delete all files from d.fielddata{idx}.Filename

    fnames_to_delete = cellfun(@getfield,diskfielddata,repmat({'Filename'},1,length(diskfielddata)),'UniformOutput',0);


    for fl = fnames_to_delete
        verbose(['Deleting ' char(fl)],30);
        delete(char(fl));
    end


else
    didx = intersect(strmatch(field,d.fieldnames),strmatch('disk',d.storagetype));
diskfielddata = {d.fielddata{didx}};
    fnames_to_delete = cellfun(@getfield,diskfielddata,repmat({'Filename'},1,length(diskfielddata)),'UniformOutput',0);


    for fl = fnames_to_delete
        verbose(['Deleting ' char(fl)],30);
        delete(char(fl));
    end

end



for fl = fnames_to_delete  %remove directory if empty
    dr = regexp(fl,'/.+/','match');

    dr = char(dr{end});
    if exist(dr,'dir')
        de = dir(dr);
        if length(de) == 2
            rmdir(dr)
        end
    end

end




    
