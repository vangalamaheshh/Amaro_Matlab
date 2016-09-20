function chunkdims = getmemchunkdims(D,field,fullfield)
%GETMEMCHUNKDIMS - Dummy function parallel to @datastruct/GETMEMCHUNKDIMS.
%
%   CHUNKDIMS = GETCHUNKDIMS(D,FIELD) Returns the size of D.(field).

%       Revisions:
%           * Function added: 10 Jan 07
%
%---
% $Id$
% $Date: 2008-02-01 17:32:51 -0500 (Fri, 01 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$



chunkdims = getsize(D,field);

    
