function s = getsize(D,field,dim)
% @DATASTRUCT/SIZE get the size of the object in FIELD of D.
%
%       S = getsize(D,FIELD,DIM)

fidx = strmatch(field,D.fieldnames,'exact');

if isempty(fidx)
    error('Field ''%s'' not found',field)
end

if strcmp('memory',D.storagetype{fidx})

    if ~exist('dim','var')
        s = size(D.fielddata{fidx});
    else
        s = size(D.fielddata{fidx},dim);
    end

elseif strcmp('disk',D.storagetype{fidx})

    if ~isempty(D.rowmapping{fidx})
        s1 = length(D.rowmapping{fidx});
    else
        s1 = D.fielddata{fidx}.Dims(1);
    end

    if ~isempty(D.colmapping{fidx})
        s2 = length(D.colmapping{fidx});
    else
        s2 = D.fielddata{fidx}.Dims(2);
    end


    if ~exist('dim','var')
        s = [s1 s2];
    elseif dim == 1
        s = s1;
    elseif dim == 2
        s = s2;
    else

        warning('Dimensions higher than 2 not currently supported for datastruct disk fields')
      s = 1;

    end
    
else
    error('Storage type not recognized')

end
