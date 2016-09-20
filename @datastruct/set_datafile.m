function M = set_datafile(M,field,fname)
%
%   SET_DATAFILE set the filename of  the hdf5 data file for a FIELD in
%   datastruct M.
%
%           M = SET_DATAFILE(M,FIELD,FILENAME)
%       
%       Revisions:
%               3 Dec 07 -- File added by Jen Dobson
%               (jdobson@broad.mit.edu)
%
%---
% $Id$
% $Date: 2008-02-28 17:02:34 -0500 (Thu, 28 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$


idx = strmatch(field,M.fieldnames,'exact');

M.fielddata{idx}.Filename = fname;

