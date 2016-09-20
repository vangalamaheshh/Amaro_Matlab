function idxout = getrowmapping(D,field,idxin)
% GETROWMAPPING returns the rowmapping of the idxin elements of field FIELD
% in D.  
%
%       IDXOUT = GETROWMAPPING(D,FIELD,IDXIN) returns the the row index for
%       that IDXIN in the associated HDF5 file.  [data is stored transposed
%       in hdf5]

%---
% $Id$
% $Rev$
% $LastChangedBy: jdobson $
% $Date: 2008-08-07 10:10:15 -0400 (Thu, 07 Aug 2008) $

fidx = strmatch(field,D.fieldnames,'exact');
idxout = D.rowmapping{fidx}(idxin);
