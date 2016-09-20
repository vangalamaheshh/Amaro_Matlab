function D = load_DWP(filename,opt)
%Temporary fix to allow a fast load of D in read-only mode

if ~exist('opt','var')
    opt = [];
end

lastwarn('');
s = load(filename,'Ds');
[a,lastid] = lastwarn;

if strcmp(lastid,'MATLAB:load:variableNotFound')
    error('Attempted to load file not saved with save_D2')
end

if strcmp(lastid,'MATLAB:unknownObjectNowStruct')
    disp('Updating class definition')
    s.Ds = datastruct(s.Ds);
end

[dum,f] = fileattrib(filename);
filename = f.Name;
fileattrib(filename,'-w');


%Find the variable that will be the loaded struct or datastruct
scell = struct2cell(s);
sclasses =  cellfun(@class,scell,'UniformOutput',0);
dsidx = strmatch('datastruct',sclasses);
cidx = strmatch('cell',sclasses);
sidx = strmatch('struct',sclasses);

if ~isempty(dsidx)
    D = scell{dsidx};
elseif ~isempty(cidx)
    D = scell{cidx};
elseif ~isempty(sidx)
    D = scell{sidx};
else
    error('No datastruct or struct found');
end

%if the datastruct class definition has changed since save of D, correct the problem!
[lastmsg,lastid] = lastwarn;
if strcmp(lastid ,'MATLAB:fieldConstructorMismatch') && isstruct(D)
    warning('Updating D to new datastruct definition')  %#ok
    XX = datastruct();
    stXX = struct(XX);
    internalfields = fieldnames(stXX);
    D = rmfield(D,setdiff(fieldnames(D),internalfields));   %take out the unneeded fields
    for fld = setdiff(internalfields,fieldnames(D))'        %add the new fields
        D.(char(fld)) = {};
    end
    D = orderfields(D,stXX);
    D = datastruct(D);   
elseif strcmp(lastid ,'MATLAB:fieldConstructorMismatch') && iscell(D)
    warning('Updating D to new datastruct definition')  %#ok
    for k = 1:length(D)
        XX = datastruct();
        stXX = struct(XX);
        internalfields = fieldnames(stXX);
        D{k} = rmfield(D{k},setdiff(fieldnames(D{k}),internalfields));  %take out the unneeded fields
        for fld = setdiff(internalfields,fieldnames(D{k}))   %add the new fields
            D{k}.(char(fld)) = {};
        end
        D{k} = orderfields(D{k},stXX);
        D{k} = datastruct(D{k});
    end
end
   



if isa(D,'cell')

    for k=1:length(D)

        if isa(D{k},'struct')
            D{k} = convert_struct(D{k});
        elseif isa(D{k},'datastruct')
            D{k} = convert_struct(D{k},opt);
            D{k} = setwriteprotect(D{k},'on');
            
            %check file pointers
            for fl = diskfieldnames(D{k})
                oldfilename = get_datafile(D{k},char(fl));
             
                if ~strcmp(oldfilename,filename)
                    D{k} = set_datafile(D{k},char(fl));
                end
            end
                            
        end

    end


elseif isa(D,'struct')

    D = convert_struct(D,opt);

elseif isa(D,'datastruct')

    D = convert_struct(D,opt);
    
    D = setwriteprotect(D,'on');

    %check file pointers
    for fl = diskfieldnames(D)
        oldfilename = get_datafile(D,char(fl));
        if ~strcmp(oldfilename,filename)    
            D =set_datafile(D,char(fl),filename);
        end

    end

end


