function state = iswriteprotected(D)
% ISWRITEPROTECTED return STATE = 1 if D is write protected; return STATE -
% 0 if D is not write protected.
%
%       STATE = ISWRITEPROTECTED(D);

%           14 Jan 07  -- jdobson@broad.mit.edu
%---
% $Id$
% $Date: 2008-07-30 14:16:42 -0400 (Wed, 30 Jul 2008) $
% $LastChangedBy: jdobson $
% $Rev$

state = getmetadata(D,'WriteProtect');
if isempty(state)
    state = 0;
end
