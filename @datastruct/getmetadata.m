function data = getmetadata(d,metaname)
%GETMETADATA get metadata from datastruct.
%
%   DATA = GETMETADATA(D,METANAME) returns metadata DATA corresponding to
%   METANAME from datastruct object D. 
%
%       See also ADDMETADATA.

%       Revisions:
%           * Function added 10 Jan 08 (jdobson@broad.mit.edu)
%---
% $Id$
% $Date: 2008-01-14 19:31:47 -0500 (Mon, 14 Jan 2008) $
% $LastChangedBy: jdobson $
% $Rev$

if ~exist('metaname','var')
    mtidx = 1:length(d.metanames);
else
mtidx = strmatch(metaname,d.metanames,'exact');
end

if isempty(mtidx)
    data = [];
else
data = d.metadata{mtidx};
end
