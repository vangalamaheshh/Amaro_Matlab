%see www.hdfgroup.com/HDF5/Tutor/selectc.html
%http://www.hdfgroup.org/HDF5/doc/index.html


%Notes:
%
%   1.  Use the user's guide reference
%   http://www.hdfgroup.org/HDF5/doc/RM_H5Front.html for interface
%   information on function calls and default settings.
%
%   2.  Need to convert to Matlab syntax:  e.g. HDFclose --> HDF.close
%
%   3.  When identifying types, use the 'NATIVE' option if possible: e.g.
%   'H5T_NATIVE_FLOAT'.  This maps to the implementation of float on our
%   system.  (which is H5T_IEEE_F32LE)
%
%   4.  Datatypes, datasets, and dataspaces must be opened and closed using
%   identifiers.  For neatness, also clear the identifiers after closing
%   their associated objects.
%
%   5.  Notice that indexing starts at 0 rather than 1, so if you want to
%   start at MATLAB's (1,1), use (0,0).  
%
%   5b.  ???  Also, C lists row size first (i.e. HOW MANY COLUMNS), then
%   column size (i.e. HOW MANY ROWS) in toplevel.Datasets.Dims.
%
%   6.  Use the HDF5info function.  Some of the lower level functions
%   require C-pointers as inputs, which matlab does not support.  So try to
%   get the information you need from HDF5 info instead.
%
%   7.  It is faster to access HDF5 rows (matlab's columns) of data than to access columns of data.
%    (Difference is a factor of 3 or 4).
%
%   8.  HUGE difference between reading data identified by coordinates and
%   data identified by continuous slabs.  Compare, for instance,
%   get_hdf5_elements with get_hdf5_cols.  Also some difference when using
%   get_hdf5_cols vs. get_hdf5_cont_cols (factor of two with size(dat) =
%   [56960 187)).
%
%   ???9.  C and MATLAB difference when indexing arrays.  For both languages,
%   reference is (col,row).  But in MATLAB, you ask, ("which elements do I
%   want as I go down this column?","which elements do I ask as I go down
%   this row?").   In C you ask ("What columns do I want?",{"What rows do I
%   want?").  


%% Ex 1: Simple load data from .mat file and write it to an HDF5 file.

cd /xchip/cancergenome/data/Rameen/jdobson/Glioma/MemIssues

load D_sim.mat

data = D.dat;


file = H5F.create('myhdf5file.h5','H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
                    %filename   %allow overwrite%

%create dataspace
dimsf = size(data);

dataspace = H5S.create_simple(ndims(D),dimsf,[]);  %last argument specifies upper limit on size of each dimension
datatype = H5T.copy('H5T_NATIVE_FLOAT');

dataset = H5D.create(file,'dataset',datatype,dataspace,'H5P_DEFAULT');
H5D.write(dataset,'H5T_NATIVE_FLOAT' ,'H5S_ALL','H5S_ALL','H5P_DEFAULT',data);
                    
     
H5D.close(dataset);
clear dataset
H5T.close(datatype);
clear datatype
H5S.close(dataspace);
clear dataspace


H5F.close(file);
clear file

fileinfo = hdf5info('myhdf5file.h5');

toplevel = fileinfo.GroupHierarchy;

%% Next: Selecting a hyperslab

%parameters required to describe hyperslab:
%       * start: starting location for the hyperslab (indexing starts at 0)
%       * stride: number of elements to separate each block to be selected.
%               [] sets the stride parameter to 1 in each dimension.
%       * count: the number of elements or blocks to select along each
%               dimension.
%       * block: the size of the blocks selected from the dataspace.  []
%               defaults to a single element in each dimension.
%
%       0 1 1 0 1 1 0 1 1 0 1 1         To select the elements with 1:
%       0 1 1 0 1 1 0 1 1 0 1 1                 start = (0,1)
%       0 1 1 0 1 1 0 1 1 0 1 1                 stride = (4,3)
%       0 0 0 0 0 0 0 0 0 0 0 0                 count = (2,4)
%       0 1 1 0 1 1 0 1 1 0 1 1                 block = (3,2)
%       0 1 1 0 1 1 0 1 1 0 1 1
%       0 1 1 0 1 1 0 1 1 0 1 1
%       0 0 0 0 0 0 0 0 0 0 0 0 


%% Ex. 2: Get a block of data from the dataset
cd /xchip/cancergenome/data/Rameen/jdobson/Glioma/MemIssues

filename = 'myhdf5file.h5';
file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT'); %open with r/w access

fileinfo = hdf5info(filename);
toplevel = fileinfo.GroupHierarchy;

datasetname = '/dataset3';
idxtl = strmatch(datasetname,{toplevel.Datasets.Name});

dataset = H5D.open(file,datasetname);

dssize = toplevel.Datasets(idxtl).Dims;

dataspace = H5D.get_space(dataset);

% Select the hyperslab of the dataset to read out
start = [0 0];  %get the 11th sample
stride = [];  %take consecutive points
count = [10 10];
block = [];
H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',start,stride,count,block);
            %this makes a hyperslab which will be used to select out from
            %the data in file

            
% Create dataspace in memory

memspace = H5S.create_simple(2,[10 10],[]); %make a column-size element in memory

% Select the hyperslab of memspace to copy in to
mstart = [0 0 ];
mstride = [];
mcount = [10 10];
mblock = [];

H5S.select_hyperslab(memspace,'H5S_SELECT_SET',mstart,mstride,mcount,mblock);

            
coldata = H5D.read(dataset,'H5T_NATIVE_FLOAT',memspace,dataspace,'H5P_DEFAULT');
%     
% H5D.close(dataset);
% clear dataset
% H5T.close(datatype);
% clear datatype
% H5S.close(dataspace);
% clear dataspace
% 
% 
% H5F.close(file);
% clear file


%% Ex. 3: Write a block of data to the dataset


newdata = single(zeros(size(coldata)));

H5D.write(dataset,'H5T_NATIVE_FLOAT',memspace,dataspace,'H5P_DEFAULT',newdata);

clear coldata

        
coldata = H5D.read(dataset,'H5T_NATIVE_FLOAT',memspace,dataspace,'H5P_DEFAULT');

%% Ex 4:  Write a 9 X 9 dataset with matlab column number as first
% numeral, matlab row number as second numeral

data2 = [1:100];
data2 = reshape(data2,10,10);
data2 = data2([1:end-1],[2:end]);

datasetname = 'dataset2';

dataspace = H5S.create_simple(2,[9 9],[]);
memspace = H5S.create_simple(2,[9 9],[]);

datasetid = H5D.create(file,datasetname,'H5T_NATIVE_INT', dataspace,'H5P_DEFAULT');

H5D.write(datasetid,'H5T_NATIVE_INT','H5S_ALL','H5S_ALL','H5P_DEFAULT',int32(data2));
H5S.close(dataspace);
H5S.close(memspace);
H5D.close(datasetid);

fileinfo = hdf5info(filename);
toplevel = fileinfo.GroupHierarchy;

datidx = strmatch(['/' datasetname],{toplevel.Datasets.Name});

%% Get the 5th column of data just written

dataset = H5D.open(file,datasetname);

dataspace = H5D.get_space(dataset);

%the dataset's hyperslab
offset = [0 4]; %get the 5th column(?)
count = [9 1];  %
H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,[],count,[]);
    %stride and block are left as []
    
%the mem dataspace
memspace = H5S.create_simple(2,[9 1],[]);

%the mem hyperslab
H5S.select_hyperslab(memspace,'H5S_SELECT_SET',[0 0],[],[9 1],[]);

coldata = H5D.read(dataset,'H5T_NATIVE_INT',memspace,dataspace,'H5P_DEFAULT');
        %% Only problem -- we got a row