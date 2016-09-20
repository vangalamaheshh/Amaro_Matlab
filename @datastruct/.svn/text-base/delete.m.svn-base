function delete(d)
%DELETE Cleanup datastruct object by deleting all hdf5 files when object is destroyed.  
%
% DELETE(D) 
disp('@datastruct/delete method')
if ~iswriteprotected(d)

    fidx = strmatch('disk',d.storagetype);

    for k = fidx'

        delete(d.fielddata{k}.Filename);

    end

end
