function d = setdirtybit(d,state)
%SETDIRTYBIT private method of datastruct.  Used to set dirtybit to ON
%during hdf5 write and turn off after write.  (Also useful because it
%touches the datastruct object on writing, indicating that a copied
%datastruct object (C = D) needs to be added to the datastruct registry)
%
% SETDIRTYBIT(D,STATE)  
%
% State can be 'on' or 'off'
%
%   jdobson@broad.mit.edu
%
%---
%$ID$
%$Date: 2008-08-07 10:09:12 -0400 (Thu, 07 Aug 2008) $
%$LastChangedBy: jdobson $
%$Rev$

if ~(strcmpi('on',state) || strcmpi('off',state))
    error('illegal state');
end

dbidx = strmatch('DirtyBit',d.metanames);
%dbdidx = strmatch(d.metanames,'DirtyBit');

if isempty(dbidx)
    d.metanames = [d.metanames 'DirtyBit'];
    d.metadata = [d.metadata 'ON'];
else
    d.metadata{dbidx} = state;
end

