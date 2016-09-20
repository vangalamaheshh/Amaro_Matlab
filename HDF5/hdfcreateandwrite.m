%sd_id = hdfsd('start','HDFmydata10.hdf','DFACC_RDWR');

sd_id = hdfsd('start','mydata10.hdf','DFACC_CREATE');

dat = D.dat;
datname = 'data2';
dattype = class(dat);
datrank = ndims(dat);
datdims = fliplr(size(dat));
%could use single and 'float' instead?
sds_id = hdfsd('create',sd_id,datname,'double',datrank,datdims);

ds_start = zeros(1:ndims(dat));
ds_stride = [];
ds_edges = fliplr(size(dat));

stat = hdfsd('writedata',sds_id,ds_start,ds_stride,ds_edges,double(dat));

stat = hdfsd('endaccess',sds_id);

stat = hdfsd('end',sd_id);


%% Read the data

sd_id = hdfsd('start','mydata10.hdf','read');

[ndatasets, nglobal_atts, stat] = hdfsd('fileinfo',sd_id);


sds_id = hdfsd('select',sd_id,1);

[dsname, dsndims, dsdims, dstype, dsatts, stat] = hdfsd('getinfo',sds_id)


ds_start = zeros(1,dsndims); % Creates the vector [0 0]
ds_stride = []; 
ds_edges = [10 10]; 

[ds_data, status] =  hdfsd('readdata',sds_id,ds_start,ds_stride,ds_edges);


stat = hdfsd('endaccess',sds_id);

stat = hdfsd('end',sd_id);


%import command:  (no need for file open and close);

%data2 = hdfread('/xchip/cancergenome/data/Rameen/jdobson/Glioma/MemIssues/mydata10.hdf', '/data2', 'Index', {[1  1],[1  1],[10  10]});