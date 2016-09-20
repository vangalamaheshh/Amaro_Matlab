function extend_hdf5_block(filename,datasetname,dat,colidx,rowidx,extend_size)
% EXTEND_HDF5_BLOCK writes data to filename, extending the dataset if
% necessary.
%
%       EXTEND_HDF5_BLOCK(FILENAME,DATASETNAME,DAT,COLIDX,ROWIDX,EXTEND_SIZ
%       E)
%
%       EXTEND_HDF5_BLOCK writes the elements of the dataset in file
%       FILENAME and dataset DATASETNAME defined by COLIDX and ROWIDX,
%       extending the dataset to EXTEND_SIZE if necessary.  DTYPE is the
%       data type of the dataset.  EXTEND_SIZE is a vector of sizes in DIM1
%       and DIM2 (matlab-style). EXTEND_HDF5_BLOCK should be used when
%       accessing a contiguous block of data or accessing partial columns
%       or rows of data.  If a contiguous block or partial rows/columns
%       cannot be found, WRITE_HDF5_ELEMENTS is called instead.
%
%

%       Revisions:
%           28 Nov 07: Function created by Jen Dobson
%           (jdobson@broad.mit.edu).
%
%           4 Dec 07:  Changed memspace selection to specify size using a
%           block (last argument) rather than count (second-to-last
%           argument) in:
%               H5S.select_hyperslab(memspace,'H5S_SELECT_SET',[0 0],[],[1 1],count);
%---
% $Id$
% $Date: 2008-08-13 13:38:31 -0400 (Wed, 13 Aug 2008) $
% $LastChangedBy: jdobson $
% $Rev$


%% Check indexing vectors

if size(colidx,1) > size(colidx,2)
    colidx = colidx';
end
if size(rowidx,1) > size(rowidx,2)
    rowidx = rowidx';
end
if size(rowidx,1) ~= 1
    error('ROWIDX must be a vector');
end
if size(colidx,1) ~= 1
    error('COLIDX must be a vector');
end

%% Sort and unique the colidx, rowidx, and dat

[colidx,csrtidx] = sort(colidx);
[rowidx,rsrtidx] = sort(rowidx);


dat = dat(rsrtidx,csrtidx);

[colidx,cunqidx,dum] = unique(colidx);
[rowidx,runqidx,dum] = unique(rowidx);

dat = dat(runqidx,cunqidx);



%% Set the data type to write to file


dtype = class(dat);

htype = mtype_to_htype(dtype);


%% Define block (select dataspace hyperslab)


%Selecting a block only works here if colidx and rowidx are sorted... at some point we
%could make this more intelligent

%If the columns are contiguous (retrieving partial rows)
if length(colidx(1):colidx(end)) == length(colidx) && isempty(find(colidx(1):colidx(end) ~= colidx))

    % First, open dataspace
    file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT'); %open with read/write access
    dataset = H5D.open(file,datasetname);
    H5D.extend(dataset,fliplr(extend_size));
    dataspace = H5D.get_space(dataset);
    

    
    %If the rows are contiguous too (i.e. the block of data is a
    %rectangle)
    if length(rowidx(1):rowidx(end)) == length(rowidx) && isempty(find(rowidx(1):rowidx(end) ~= rowidx))

        offset = [colidx(1)-1 rowidx(1)-1];
        stride = [];
        count =  [1 1];
        block =  [colidx(end)-colidx(1)+1  rowidx(end)-rowidx(1)+1];        %dimensions of block of data

        H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,stride,count,block);


        %If the rows aren't also contiguous (i.e. get partial rows) or if the
        %rowidx aren't sorted
    else

        offset = [colidx(1)-1 rowidx(1)-1];
        stride = [];
        count =  [1 1];
        block =  [colidx(end)-colidx(1)+1 1];        %block of one row

        H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,stride,count,block);
        %append to the dataset's hyperslab
        for k = rowidx(2:end)
            offset = [colidx(1)-1 k-1];
            H5S.select_hyperslab(dataspace,'H5S_SELECT_OR',offset,stride,count,block);
        end

    end

    %make memory space 
    count = size(dat);
    %the mem dataspace
    memspace = H5S.create_simple(2,count,[]);
    %the mem hyperslab
    H5S.select_hyperslab(memspace,'H5S_SELECT_SET',[0 0],[],[1 1],count);
    %% Write data
    H5D.write(dataset,htype,memspace,dataspace,'H5P_DEFAULT',dat);
    %% Close resources
    H5S.close(dataspace);
    H5D.close(dataset);
    H5F.close(file)

    %If the rows are contiguous but not the cols (i.e. get partial cols) or if the colidx aren't sorted
elseif length(rowidx(1):rowidx(end)) == length(rowidx) && isempty(find(rowidx(1):rowidx(end) ~= rowidx))
    
    %open dataspace
    file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT'); %open with read only access
    dataset = H5D.open(file,datasetname);
    H5Dextend(dataset,fliplr(extend_size));
    dataspace = H5D.get_space(dataset);

    offset = [colidx(1)-1 rowidx(1)-1];
    stride = [];
    count =  [1 1];
    block =  [1 rowidx(end)-rowidx(1)+1];        %block of one row

    H5S.select_hyperslab(dataspace,'H5S_SELECT_SET',offset,stride,count,block);
    %append to the dataset's hyperslab
    for k = colidx(2:end)
        offset = [k-1 rowidx(1)-1];
        H5S.select_hyperslab(dataspace,'H5S_SELECT_OR',offset,stride,count,block);
    end

    %make memory space
    count = size(dat);
    %the mem dataspace
    memspace = H5S.create_simple(2,count,[]);
    %the mem hyperslab
    H5S.select_hyperslab(memspace,'H5S_SELECT_SET',[0 0],[],[1 1],count);
    %% Write data
    H5D.write(dataset,htype,memspace,dataspace,'H5P_DEFAULT',dat);
    %% Close resources
    H5S.close(dataspace);
    H5D.close(dataset);
    H5F.close(file)





    %%%%%%If can't make a block, then do a recursive write on smaller
    %%%%%%blocks
else
    warning('Could not make single block of referenced indices.  Writing multiple blocks.')

    [rowbounds,dridx] = findboundaries(rowidx);
    [colbounds,dcidx] = findboundaries(colidx);

    for j = 1:length(rowbounds)
        for k = 1:length(colbounds)

            extend_hdf5_block(filename,datasetname,dat(dridx{j}(1):dridx{j}(2),dcidx{k}(1):dcidx{k}(2)),...
                colbounds{k}(1):colbounds{k}(2),rowbounds{j}(1):rowbounds{j}(2),extend_size);

        end
    end


end






end

function [B,idx] = findboundaries(vect);
%vector v must be a sorted vector

cmpset = setdiff(min(vect):max(vect),vect);
v = vect;
k = 1;
B = {};
idx = {};
while ~isempty(v)

    B{k}(1) = v(1);  %v's minimum
    v = v(2:end);  %strip off minimum
    if ~isempty(cmpset)
        Bk2 = max(v(find(v<min(cmpset))));
    else
        Bk2 = max(v);
    end

    if isempty(Bk2)
        B{k}(2) = B{k}(1);
    else
        B{k}(2) = Bk2;
    end

    v = v(v>B{k}(2));

    if ~isempty(v)
        cmpset = cmpset(cmpset>min(v));
    else
        cmpset = [];
    end


    idx{k}(1) = find(vect == B{k}(1));
    idx{k}(2) = find(vect == B{k}(2));
    k = k+1;


end

end




    
    
    
    

