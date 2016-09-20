function filename = get_datafile(M,field)
%
%   GET_DATAFILE get the filename of  the hdf5 data file for a FIELD in
%   datastruct M.
%
%           FILENAME = GET_DATAFILE(M,FIELD)
%       
%       Revisions:
%               3 Dec 07 -- File added by Jen Dobson
%               (jdobson@broad.mit.edu)
%
%---
% $Id$
% $Date: 2008-02-14 09:17:44 -0500 (Thu, 14 Feb 2008) $
% $LastChangedBy: jdobson $
% $Rev$


idx = strmatch(field,M.fieldnames,'exact');

filename = M.fielddata{idx}.Filename;


