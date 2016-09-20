function rewrite_dset_attribute(file,datasetname,attrname,val)
%REWRITE_DSET_ATTRIBUTE write an attribute to data set.
%   
%   REWRITE_DSET_ATTRIBUTE(FILE,DATASETNAME,ATTRNAME,VAL)

%---
% $Id$
% $Date: 2007-11-30 14:17:30 -0500 (Fri, 30 Nov 2007) $
% $LastChangedBy: jdobson $
% $Rev$


%don't use full name
pat = '/.*/';
attrname = regexprep(attrname,pat,'');

fid = H5F.open(file,'H5F_ACC_RDWR','H5P_DEFAULT');
dsid = H5D.open(fid,datasetname);

attrid = H5A.open_name(dsid,attrname);

H5A.write(attrid,'H5ML_DEFAULT',val);

H5A.close(attrid);
H5D.close(dsid);
H5F.close(fid);