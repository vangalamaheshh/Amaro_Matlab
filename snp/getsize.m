function s = getsize(D,field,dim)
%SIZE get the size of the object in FIELD of D.
%
%       S = getsize(D,FIELD,DIM)

if ~exist('dim','var')    
    s = size(D.(field));
else
    s = size(D.(field),dim);
end
