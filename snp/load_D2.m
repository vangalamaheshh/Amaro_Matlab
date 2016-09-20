function D=load_D2(fname,hdf5dir,varargin)
% LOAD_D2 load a datastruct object.
%
%   D = LOAD_D2(FNAME) loads a datastruct object from the
%   file FNAME.  FNAME is a .mat file created by SAVE_D2.  
%
%   D = LOAD_D2(FNAME,HDF5DIR)  When the object
%   is loaded, the hdf5 files associated with the diskfields from that
%   object will be created in HDF5DIR.  If HDF5DIR is not specified, pwd
%   will be used. 
%
%   D = LOAD_D2(FNAME,HDF5DIR,'skip-cell') To skip the conversion of cell-of-cell fields (such as
%   used-normals), set the third argument to 'skip-cell'.  HDF5DIR may
%   optionally be left blank.
%
%       See also: save_D2

%       Revisions:
%            22 Dec 07 -- File created (jdobson@broad.mit.edu)
%
%---
% $Id$
% $Date: 2008-08-12 12:14:38 -0400 (Tue, 12 Aug 2008) $
% $LastChangedBy: jdobson $
% $Rev$

if ~isempty(varargin) && ischar(varargin{1})
    if strcmp(varargin{1},'skip-cell')
        opt = varargin{1};
    end
else
    opt = [];
end


lastwarn('');
load(fname,'Ds');
[lastmsg,lastid] = lastwarn;

if strcmp(lastid,'MATLAB:load:variableNotFound')
    error('Attempted to load file not saved with save_D2')
end

if strcmp(lastid,'MATLAB:unknownObjectNowStruct')
    disp('Updating class definition')
    Ds = datastruct(Ds);
end


D = Ds;
clear Ds;

% %Find the variable that will be the loaded struct or datastruct
% scell = struct2cell(s);
% sclasses =  cellfun(@class,scell,'UniformOutput',0);
% dsidx = strmatch('datastruct',sclasses);
% cidx = strmatch('cell',sclasses);
% sidx = strmatch('struct',sclasses);
% 
% if ~isempty(dsidx)
%     D = scell{dsidx};
% elseif ~isempty(cidx)
%     D = scell{cidx};
% elseif ~isempty(sidx)
%     D = scell{sidx};
% else
%     error('No datastruct or struct found');
% end

%if the datastruct class definition has changed since save of D, correct the problem!

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
    warning('Updating D to new datastruct definition') %#ok
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
   
    
% Set the hdf5 directory
if ~exist('hdf5dir','var') || isempty(hdf5dir)
    hdf5dir = [pwd filesep];
elseif ~strcmp('file_sep',hdf5dir(end))
    hdf5dir = [hdf5dir filesep];
end


% if length(fieldnames(s))>1
%     error('file should contain only one variable')
% end

%set filename to full filename
[s,f] = fileattrib(fname);
fname = f.Name;


if isa(D,'cell')

    for k=1:length(D)

        if isa(D{k},'struct')
            D{k} = convert_struct(D{k});
        elseif isa(D{k},'datastruct')
            diskfields = diskfieldnames(D{k})';


            D{k} = convert_struct(D{k},opt);

            if ~isempty(diskfields)
                for fl = setdiff(fieldnames(D{k}),diskfields)'
                    if getsize(D{k},char(fl)) == getsize(D{k},diskfields{1})
                        disp(['Converting field ''' char(fl) ''' to diskfield'])
                        D{k} = convert_to_diskfield(D{k},char(fl),get_new_filenames(strcat(hdf5dir,char(fl))),char(fl));         
                    end
                end
            end

      
           
            D{k} = change_data_pointer(D{k},0,diskfields,strcat(hdf5dir,diskfields,'_',num2str(k),'.h5'),fname);

        end

    end


elseif isa(D,'struct')

    D = convert_struct(D,opt);

elseif isa(D,'datastruct')
    diskfields = diskfieldnames(D)';

    D = convert_struct(D,opt);
    if ~isempty(diskfields)
        for fl = setdiff(fieldnames(D),diskfields)'
            if isequal(getsize(D,fl), getsize(D,diskfields(1)))
                D = convert_to_diskfield(D,char(fl),get_new_filenames(strcat(hdf5dir,char(fl),'.h5')),char(fl));
            
            end
        end
    end
  
    D = change_data_pointer(D,0,diskfields,strcat(hdf5dir,diskfields,'.h5'),fname);

end
end

