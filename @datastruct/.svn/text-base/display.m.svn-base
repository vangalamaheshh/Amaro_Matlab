function display(D)
%D/DISPLAY  Command Window Display of a D.

disp(' ');
disp([inputname(1),' = '])
disp(' ');

%
% disp(D.memfields);
%
% if ~isempty(D.diskfields)
% diskfieldnames = fieldnames(D.diskfields);

% If array, describe contents of array
if ~isequal(size(D), [1 1])
    fprintf(1,'%dx%d datastruct array\n\n',size(D,1),size(D,2))
    return
end

% If empty, show empty
if isempty(D.fielddata)
    fprintf(1,'%7s\n','[]');
    return
end
   
memfldidx = strmatch('memory',D.storagetype);

for k = memfldidx'
    
    if isnumeric(D.fielddata{k})

        str = ['%14s: [%d' repmat('x%d',1,ndims(D.fielddata{k})-1) ' %s]\n'];
        fprintf(1,str,D.fieldnames{k},size(D.fielddata{k}),class(D.fielddata{k}))

    elseif ischar(D.fielddata{k})
        if ndims(D.fielddata{k} == 2) && size(D.fielddata{k},1)== 1 ...
                && size(D.fielddata{k},2) < 100

            fprintf(1,'%14s: ''%s''\n',D.fieldnames{k},D.fielddata{k})

        else

            str = ['%14s: [%d' repmat('x%d',1,ndims(D.fielddata{k})-1) ' %s]\n'];
            fprintf(1,str,D.fieldnames{k},size(D.fielddata{k}),class(D.fielddata{k}))
        end

    elseif isstruct(D.fielddata{k})
        fprintf(1,'%14s: [%dx%d %s]\n',D.fieldnames{k},size(D.fielddata{k},1),size(D.fielddata{k},2),class(D.fielddata{k}))

    elseif  iscell(D.fielddata{k})
        fprintf(1,'%14s: {%dx%d %s}\n',D.fieldnames{k},size(D.fielddata{k},1),size(D.fielddata{k},2),class(D.fielddata{k}))

    else
        warning('Unknown fieldtype in Datastruct')

    end
end

datfldidx = strmatch('disk',D.storagetype);



for k=datfldidx'
    
    if isfield(D.fielddata{k},'Dims') && length(D.fielddata{k}.Dims) == 2
        
        if ~isempty(D.colmapping{k})
            d2 = length(D.colmapping{k});
        else
            d2 = D.fielddata{k}.Dims(2);
        end
        
        if ~isempty(D.rowmapping{k})
            d1 = length(D.rowmapping{k});
        else
            d1 = D.fielddata{k}.Dims(1);
        end
                   
        
        
        fprintf(1,'%14s: [%dx%d %s]\n',D.fieldnames{k},d1,d2,D.fielddata{k}.Datatype.Class)
        
      
        
        
    elseif isfield(D.fielddata{k},'Dims') && length(D.fielddata{k}.Dims == 1)
        
        if ~isempty(D.colmapping{k})
            d1 = length(D.colmapping{k});
        elseif ~isempty(D.rowmapping{k})
            d1 = length(D.rowmapping{k});
        else
            d1 = D.fielddata{k}.Dims(1);
       
        end
        
            
        
        fprintf(1,'%14s: [%dx1 %s]\n',D.fieldnames{k},d1,D.fielddata{k}.Datatype.Class)
        
    else
        
        fprintf(1,'%14s: ???\n',D.fieldnames{k})
    end

end

