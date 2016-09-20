function merged = mergestruct(varargin)
%MERGESTRUCTCombine structures.
% 
%   merged = mergestruct(varargin)
%
% collect field names
flds = [];
 for k = 1:nargin
    try
        flds = [flds ; fieldnames(varargin{k})];
    catch MEstruct
        throw(MEstruct)
    end
end

% ensure the field names are unique.
if length(flds) ~= length(unique(flds))
    error('mergestruct:FieldsNotUnique',...
        'Field names must be unique');
end
    
% concatenate the data from each struct.
c = [];
for k = 1:nargin
    try
        c = [c ; struct2cell(varargin{k})];
    catch MEdata
        throw(MEdata);
    end
end
% Construct the output.
merged = cell2struct(c, flds, 1);
