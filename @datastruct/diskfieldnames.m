function nms = diskfieldnames(D)
%   DISKFIELDNAMES list the diskfields of datastruct D and return in cell
%   array NMS.
%
%       

%           Revisions:
%               29 Nov 07 - Added (jdobson@broad.mit.edu)
%---
% $Id$
% $Date: 2007-11-30 13:47:21 -0500 (Fri, 30 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$



stor = D.storagetype;
diskidx = strmatch('disk',stor,'exact');

nms = D.fieldnames(diskidx);
