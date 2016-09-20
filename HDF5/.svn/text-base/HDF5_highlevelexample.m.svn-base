

%% HDF5 High Level Example
cwd = pwd;

cd /xchip/cancergenome/data/Rameen/jdobson/Glioma/MemIssues

load D_sim.mat

cd(pwd)

dat = D.dat;

whos dat
%   Name           Size                Bytes  Class     Attributes
%   dat       249989x100            99995600  single  

hdf5write('myhdf5file.h5','dataset1',dat)

fileinfo = hdf5info('myhdf5file.h5')

% fileinfo = 
%           Filename: 'myhdf5file.h5'
%         LibVersion: '1.6.5'
%             Offset: 0
%           FileSize: 100000264
%     GroupHierarchy: [1x1 struct]

toplevel = fileinfo.GroupHierarchy

% toplevel = 
%       Filename: 'myhdf5file.h5'
%           Name: '/'
%         Groups: []
%       Datasets: [1x1 struct]
%      Datatypes: []
%          Links: []
%     Attributes: []

toplevel.Datasets

% ans = 
%       Filename: 'myhdf5file.h5'
%           Name: '/dataset1'
%           Rank: 2
%       Datatype: [1x1 struct]
%           Dims: [249989 100]
%        MaxDims: [249989 100]
%         Layout: 'chunked'
%     Attributes: []
%          Links: []
%      Chunksize: []
%      FillValue: 0

clear

dat = hdf5read('myhdf5file.h5','/dataset1');

whos
%   Name           Size                Bytes  Class     Attributes
%   dat       249989x100            99995600  single        
