function [fn,ds,fieldnames] = gethdf5info_lite(D)
%GETHDF5INFO_LITE returns the hdf5 filename and dataset name for all disk
%fields in D.
%
%       [FN,DS,FIELDNAMES] = GETHDF5INFO_LITE(D)
%
%   added 26 FEB 08  (jdobson)

%--------
%$Id$
%$Date: 2008-08-13 13:35:27 -0400 (Wed, 13 Aug 2008) $
%$LastChangedBy: jdobson $
%$Rev$

didx = strmatch('disk',D.storagetype);
fielddata = [D.fielddata{didx}];
if ~isempty(fielddata)
fn = {fielddata.Filename};
ds = {fielddata.Name};

fieldnames = D.fieldnames(didx);
else
    fn = [];
    ds = [];
    fieldnames = [];
end

