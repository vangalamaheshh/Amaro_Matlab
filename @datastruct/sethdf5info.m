function D = sethdf5info(D,fields,hdf5field,setto)
%SETHDF5INFO sets the hdf5 info in a datastruct.  This function does not
%operate on the HDF5 file itself!!  (Use with care or you could destroy your link to the file.)
%
%   D = SETHDF5INFO(D,FIELDS,HDF5FIELD,SETTO)
%
%   added 26 FEB 08  (jdobson)

%--------
%$Id$
%$Date: 2008-02-29 12:45:32 -0500 (Fri, 29 Feb 2008) $
%$LastChangedBy: jdobson $
%$Rev$
if ~iscell(setto)
    setto = {setto};
end


[iidx,didx,fidx]=intersect(D.fieldnames,fields);
setto = setto(fidx);

fielddata = D.fielddata(didx);

fdnames = cellfun(@fieldnames,fielddata,'UniformOutput',0);
hdf5idx = cellfun(@strmatch,repmat({hdf5field},1,length(fdnames)),fdnames);
fielddata = cellfun(@struct2cell,fielddata,'UniformOutput',0);

fielddata = horzcatfill(fielddata{:});
fielddata = [fielddata{:}];
offset = repmat(size(fielddata,1),1,size(fielddata,2));
offset = cumsum(offset);
offset = [0 offset(1:(end-1))];
hdf5idx = offset + hdf5idx;
fielddata(hdf5idx) = setto;

%convert fielddata from nXm cell to nX1 cell
fun = @(x) fielddata(:,x);
vect = mat2cell(1:size(fielddata,2),1,ones(1,size(fielddata,2)));
fielddata = cellfun(fun,vect,'UniformOutput',0);

fielddata = cellfun(@cell2struct,fielddata,fdnames,'UniformOutput',0);


D.fielddata(didx) = fielddata;