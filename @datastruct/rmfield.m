function s = rmfield(s,field)
%RMFIELD Remove fields from a datastruct array.
%   S = RMFIELD(S,field) removes the specified field from the
%   datastruct array S. 
%
%   S = RMFIELD(S,FIELDS) removes more than one field at a time
%   when FIELDS is a character array or cell array of strings.  The
%   changed structure is returned. The size of input S is preserved.
%

%

%--------------------------------------------------------------------------------------------
% handle input arguments
if ~isa(s,'datastruct') 
    error('S must be a DATASTRUCT object.'); 
end

if ~ischar(field) && ~iscellstr(field)
   error( 'FIELDNAMES must be a string or a cell array of strings.');
elseif ischar(field)
   field = cellstr(field); 
end

% get fieldnames of struct
f = s.fieldnames;

% Determine which fieldnames to delete.
idxremove = [];
for i=1:length(field)
  j = strmatch(field{i},f,'exact');
  if isempty(j)
    if length(field{i}) > namelengthmax
      error('The string ''%s''\nis longer than the maximum MATLAB name length.  It is not a valid fieldname.',...
          field{i});
    else
      error('MATLAB:rmfield:InvalidFieldname', 'A field named ''%s'' doesn''t exist.',field{i});
    end
  end
  idxremove = [idxremove;j];
end

% set indices of fields to keep
idxkeep = 1:length(f);
idxkeep(idxremove) = [];

if ~iswriteprotected(s)
     % delete hdf5 files associated with fields to remove
    idxdelete = intersect(idxremove,strmatch('disk',s.storagetype));
    if ~isempty(idxdelete)
        for k = idxdelete
            delete(s.fielddata{k}.Filename);
        end
    end
end

for stfld = {'fieldnames','fielddata','storagetype','colmapping','rowmapping'}
    s.(char(stfld)) = s.(char(stfld))(idxkeep);
end



%--------------------------------------------------------------------------------------------
