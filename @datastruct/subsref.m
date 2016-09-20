function b=subsref(a,s)
%@DATASTRUCT/SUBSREF overloaded method to reference contents of class DATASTRUCT
%
%       20 Nov 07: Created by Jen Dobson (jdobson@broad.mit.edu)
%
%       27 Nov 07:  Updated to read hdf5 files using blocked data if
%       possible.

%---
% $Id$
% $Date: 2011-04-11 15:22:11 -0400 (Mon, 11 Apr 2011) $
% $LastChangedBy: schum $
% $Rev$


% If reference is to a field of the datastruct
if strmatch(s(1).type,'.')

    idx = strmatch(s(1).subs,a.fieldnames,'exact');

    if isempty(idx)
        error('Field name does not exist')
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Memory field  (easy case)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strmatch(a.storagetype(idx),'memory')

        if length(s) == 1 && strcmp('.',s(1).type)
            b = a.fielddata{idx};

        elseif length(s) >= 2 && strcmp('.',s(1).type)
            if isstruct(a.fielddata{idx}) && length(s) == 2 && s(2).type(1) == '.'


                b  = {a.fielddata{idx}.(s(2).subs)};


            else
                b = subsref(a.fielddata{idx},s(2:end));
            end
        else
            error('Illegal reference')
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Disk field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strmatch(a.storagetype(idx),'disk')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% READ DATA SUBSET IF LENGTH OF S is 2, then need subset
        %%%%%%%%%%% of the field. (indicated by subscripted references,
        %%%%%%%%%%% logicals, logicals and subscripted, or linear indices)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% If reading data subset, and s(2) has length 2 (==> x,y
        %%%%%%%%%%% referenced independently with subscripted references or
        %%%%%%%%%%% logicals)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(s) == 2 && length(s(2).subs)==2 && strcmp('()',s(2).type) && strcmp('.',s(1).type)
            cols = s(2).subs{2};
            rows = s(2).subs{1};

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% PARSE COLUMN INPUT (whether logical, subsripted, or
            %%%%%% colon
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isnumeric(cols)

                if isempty(cols)
                    b = [];
                    return
                end
                
                if min(cols) < 1 || max(cols) > getsize(a,a.fieldnames{idx},2)
                    error('Index exceeds dimensions')
                end


                if size(cols,1) > 1
                    if size(cols,2) == 1
                        cols = cols';
                    else
                        error('Indices must be vectors')
                    end
                end

            elseif ischar(cols) && strcmp(cols,':')

                cols = 1:getsize(a,a.fieldnames{idx},2);

            elseif islogical(cols)

                cols = find(cols);
                
                if isempty(cols)
                    b = [];
                    return
                end
                
                if min(cols) < 1 || max(cols) > getsize(a,a.fieldnames{idx},2)
                    error('Index exceeds dimensions')
                end


                if size(cols,1) > 1
                    if size(cols,2) == 1
                        cols = cols';
                    else
                        error('Indices must be vectors')
                    end
                end

            else

                error('Unrecognized subscript for column indices')

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% PARSE ROW INPUT (whether logical, subsripted, or
            %%%%%% colon
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if isnumeric(rows)

                if isempty(rows)
                    b = [];
                    return
                end
                    
                if min(rows) < 1 || max(rows) > getsize(a,a.fieldnames{idx},1)
                    error('Index exceeds dimensions')
                end

                if size(rows,1) > 1
                    if size(rows,2) == 1
                        rows = rows';
                    else
                        error('Indices must be vectors')
                    end
                end


            elseif ischar(rows) && strcmp(rows,':') && ~ischar(cols)

                rows = 1:getsize(a,a.fieldnames{idx},1);

            elseif islogical(rows)

                rows = find(rows);

                if isempty(rows)
                    b = [];
                    return
                end
                
                if min(rows)<1 || max(rows) > getsize(a,a.fieldnames{idx},1)
                    error('Index exceeds dimensions')
                end

                if size(rows,1) > 1
                    if size(rows,2) == 1
                        rows = rows';
                    else
                        error('Indices must be vectors')
                    end
                end

            else

                error('Unrecognized subscript for row indices')

            end



            %verbose('Reading Data ''%s''from File ''%s''\n',30,a.fielddata{idx}.Name,a.fielddata{idx}.Filename);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % map to colmapping and rowmapping
            %%%%%%%%
            %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(a.colmapping{idx})
                dc = a.colmapping{idx}(cols);
            else
                dc = cols;
            end
            if ~isempty(a.rowmapping{idx})
                dr = a.rowmapping{idx}(rows);
            else
                dr = rows;
            end
    

            if isequal(cols,1:getsize(a,a.fieldnames{idx},2)) && isequal(rows,getsize(a,a.fieldnames{idx},1))

                b = hdf5read(a.fielddata{idx}.Filename,a.fielddata{idx}.Name);
                b = b(dr,dc);
            else

                b = get_hdf5_block(a.fielddata{idx}.Filename,a.fielddata{idx}.Name,dc,dr,a.fielddata{idx}.Datatype.Class);
            end
            
     
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% If reading data subset, and s(2) has length 1
        %%%%%%%%%%% (i.e. logical or linear index referencing)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        elseif length(s) == 2 && length(s(2).subs) == 1 && strcmp('()',s(2).type) && strcmp('.',s(1).type) ...
                && getsize(a,a.fieldnames{idx},1) > 1 && getsize(a,a.fieldnames{idx},2) > 1    %%% Last 2 check that the field is 2-dim array

            %If logical array is passed
            if islogical(s(2).subs{1})

                if size(s(2).subs{:},1) ~= getsize(a,a.fieldnames{idx},1) || size(s(2).subs{:},2) ~= getsize(a,a.fieldnames{idx},2)
                    error('Input logical reference array has dimensions different from value array')
                else
                    [rowcoords,colcoords] = ind2sub(getsize(a,a.fieldnames{idx}),find(s(2).subs{:}));
                end

            %If linear index to 2 dimensional data is passed
            else

                [rowcoords,colcoords] = ind2sub(getsize(a,a.fieldnames{idx}),s(2).subs{:});

            end

            %verbose('Reading Data ''%s''from File ''%s''\n',30,a.fielddata{idx}.Name,a.fielddata{idx}.Filename);

            if ~isempty(a.colmapping{idx})
            dc = a.colmapping{idx}(colcoords);
            else
                dc = colcoords;
            end
            
            if ~isempty(a.colmapping{idx})
            dr = a.rowmapping{idx}(rowcoords);
            else
                dr = rowcoords;
            end
            


            b = get_hdf5_elements(a.fielddata{idx}.Filename,a.fielddata{idx}.Name,dc,dr,a.fielddata{idx}.Datatype.Class);

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%  If 1-dimensional array
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif length(s) == 2 && length(s(2).subs) == 1 && strcmp('()',s(2).type) && strcmp('.',s(1).type) ...
                && (getsize(a,a.fieldnames{idx},1) == 1 || getsize(a,a.fieldnames{idx},2) == 1 )    %%% Last 2 check that the field is col or row vector

            %If logical array is passed
            if islogical(s(2).subs{1})

                if size(s(2).subs{:}) == getsize(a,a.fieldnames{idx},1) 
                    rowcoords = find(s(2).subs{1});
                    colcoords = ones(1,length(rowcoords));
                elseif size(s(2).subs{:}) == getsize(a,a.fieldnames{idx},2)
                    colcoords = find(s(2).subs{1});
                    rowcoords = ones(1,length(colcoords));
                else
                    error('Input logical reference array has dimensions different from value array')
                end

                %If linear index or ':' is passed
            elseif isnumeric(s(2).subs{1})
                
                if getsize(a,a.fieldnames{idx},1) == 1
                    colcoords = s(2).subs{1};
                    rowcoords = ones(1,length(colcoords));
                elseif getsize(a,a.fieldnames{idx},2) == 1
                    rowcoords = s(2).subs{1};
                    colcoords = ones(1,length(rowcoords));
                end
                
            elseif strcmp(':',s(2).subs{1})
                
                if getsize(a,a.fieldnames{idx},1) == 1
                    colcoords = 1:getsize(a,a.fieldnames{idx},2);
                    rowcoords = ones(1,length(colcoords));
                elseif getsize(a,a.fieldnames{idx},2) == 1
                    rowcoords = 1:getsize(a,a.fieldnames{idx},1);
                    colcoords = ones(1,length(rowcoords));
                end

            end
            
            if ~isempty(a.colmapping{idx})
                dc = a.colmapping{idx}(colcoords);
            else
                dc = colcoords;
            end

            if ~isempty(a.colmapping{idx})
                dr = a.rowmapping{idx}(rowcoords);
            else
                dr = rowcoords;
            end



            b = get_hdf5_elements(a.fielddata{idx}.Filename,a.fielddata{idx}.Name,dc,dr,a.fielddata{idx}.Datatype.Class);


            %verbose('Reading Data ''%s''from File ''%s''\n',30,a.fielddata{idx}.Name,a.fielddata{idx}.Filename);

 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%  NO SUBSCRIPTS: Read all data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif length(s) == 1 && strcmp('.',s(1).type)
    
            %verbose('Reading Data ''%s''from File ''%s''\n',30,a.fielddata{idx}.Name,a.fielddata{idx}.Filename);
            b = hdf5read(a.fielddata{idx}.Filename,a.fielddata{idx}.Name);
            
            if ~isempty(a.colmapping{idx})
            dc = a.colmapping{idx};
            else
                dc = 1:a.fielddata{idx}.Dims(2);
            end
            
            if ~isempty(a.rowmapping{idx})
                dr = a.rowmapping{idx};
            else
                dr = 1:a.fielddata{idx}.Dims(1);
            end
            
             
            
            b = b(dr,dc);

            
        else
            error('Illegal reference.')
        end
    else
        % a.storagetype(idx) neither memory nor disk
        error('Unrecognized storage type')
    end

    
    % If reference is to an index in an array of datastructs
elseif strmatch(s(1).type,'()')
    if length(s) == 1
        b = a(s(1).subs{:});
    else
        b = subsref(a(s(1).subs{:}),s(2:end));
    end
    
else
    error('Datastruct object must be referenced with . or ()')
end


