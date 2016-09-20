function D = setwriteprotect(D,state)
%   SETWRITEPROTECT -- toggle or set the writeprotect state of D.
%
%       D = SETWRITEPROTECT(D,STATE)
%
%           STATE is 'ON' or 'OFF'

%---
% $Id$
% $Date: 2008-02-15 14:28:00 -0500 (Fri, 15 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$


if strcmpi('ON',state)
    D = setmetadata(D,'WriteProtect',1);
elseif strcmpi('OFF',state)
    D = setmetadata(D,'WriteProtect',0);
end
