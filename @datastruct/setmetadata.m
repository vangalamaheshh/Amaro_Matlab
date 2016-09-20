function d = setmetadata(d,metaname,data)
%SETMETADATA set metadata in datastruct.
%
%   D = SETMETADATA(D,METANAME,DATA) adds the metadata DATA specified by string
%   METANAME to the datastruct object D.  Possible fields are: MemLimit
%   (max number of bytes to store in memory on data transfers), DirtyBit,
%   ReadOnly.

%       Revisions:
%           * Function added 10 Jan 07 (jdobson@broad.mit.edu)
%
%---
% $Id$
% $Date: 2008-01-14 19:26:42 -0500 (Mon, 14 Jan 2008) $
% $LastChangedBy: jdobson $
% $Rev$


idx = strmatch(metaname,d.metanames,'exact');

if isempty(idx)
    idx = length(d.metanames) + 1;
end

d.metadata{idx} = data;
d.metanames{idx} = metaname;
