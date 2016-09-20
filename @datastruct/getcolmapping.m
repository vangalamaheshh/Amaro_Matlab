function idxout = getcolmapping(D,field,idxin)
% GETCOLMAPPING returns the colmapping of the idxin elements of field FIELD
% in D.  
%
%       IDXOUT = GETCOLMAPPING(D,FIELD,IDXIN) returns the the row index for
%       that IDXIN in the associated HDF5 file.  [data is stored transposed
%       in hdf5]

%---
% $Id$
% $Rev$
% $LastChangedBy: jdobson $
% $Date: 2008-02-28 17:08:45 -0500 (Thu, 28 Feb 2008) $

fidx = strmatch(field,D.fieldnames,'exact');
idxout = D.colmapping{fidx}(idxin);
