function [D,datmean] = remove_nans(D,rthresh,cthresh)
%REMOVE_NANS remove nans from unsegmented copy number data
%
%   [D,DATMEAN] = remove_nans(D,RTHRESH,CTHRESH)
%
% REMOVE_NANS removes all NaNs from D.dat. If a marker has a fraction of
% samples that are NaNs exceeding RTHRESH, the marker will be removed from
% D. If a sample has greater than the fraction CTHRESH of its markers as
% NaNs, the sample will be removed from D. Any remaining NaNs are replaced
% with the calculated value 
%
%              Row.average * Col.average / DATMEAN
% 
% Where Row.average and Col.average are nanmeans across samples for each
% marker and markers for each sample respectively, and DATMEAN is the
% overall mean value (optional second return value).

%% process datastructs passed in a cell array one at a time
if iscell(D)
    datmean = nan(size(D)); %! or weighted mean?
    for i=1:numel(D)
        [D{i} datmean(i)] = remove_nans(D{i},rthresh,cthresh);
    end
    return
end

%% eliminate rows with too many NaNs
datsize = getsize(D,'dat');
nans_in_row = zeros(datsize(1),1);

chunkdims = getmemchunkdims(D,'dat');
% loop over row chunks
for i = 1:chunkdims(1):datsize(1)
    % create row chunk-selection index
    idx1 = i:min(i+chunkdims(1)-1,datsize(1));
    % loop over column chunks
    for j = 1:chunkdims(2):datsize(2)
        % construct logical column index
        idx2 = j:min(j+chunkdims(2)-1,datsize(2));
        % get marginal count of nans in the chunk
        dat = D.dat(idx1,idx2);
        nans_in_row(idx1) = nans_in_row(idx1) + sum(isnan(dat),2);
    end
end
D = reorder_D_rows(D,nans_in_row/datsize(2) <= rthresh);

%% eliminate columns with too many NaNs
datsize = getsize(D,'dat');
nans_in_col = zeros(1,datsize(2));

chunkdims = getmemchunkdims(D,'dat');
% loop over row chunks
for i = 1:chunkdims(1):datsize(1)
    % create row chunk-selection index
    idx1 = i:min(i+chunkdims(1)-1,datsize(1));
    % loop over column chunks
    for j = 1:chunkdims(2):datsize(2)
        % construct logical column index
        idx2 = j:min(j+chunkdims(2)-1,datsize(2));
        % get marginal count of nans in the chunk
        dat = D.dat(idx1,idx2);
        nans_in_col(idx2) = nans_in_col(idx2) + sum(isnan(dat),1);
    end
end
D = reorder_D_cols(D,nans_in_col/datsize(1) <= cthresh);

%% calculate marginal nanmeans 
datsize = getsize(D,'dat');

rowsums = zeros(datsize(1),1);
colsums = zeros(1,datsize(2));
nums_in_row = zeros(datsize(1),1);
nums_in_col = zeros(1,datsize(2));

% loop over row chunks
for i = 1:chunkdims(1):datsize(1)
    % create row chunk-selection index
    idx1 = i:min(i+chunkdims(1)-1,datsize(1));
    % loop over column chunks
    for j = 1:chunkdims(2):datsize(2)
        % construct column chunk-selection index
        idx2 = j:min(j+chunkdims(2)-1,datsize(2));
        % get chunk
        dat = D.dat(idx1,idx2);
        % calculate marginal nansums for chunk
        rowsums(idx1) = rowsums(idx1) + nansum(dat,2);
        colsums(idx2) = colsums(idx2) + nansum(dat,1);
        % count non-nans in chunk
        nums_in_row(idx1) = nums_in_row(idx1) + sum(~isnan(dat),2);
        nums_in_col(idx2) = nums_in_col(idx2) + sum(~isnan(dat),1);
    end
end
rowmeans = rowsums ./ nums_in_row;
colmeans = colsums ./ nums_in_col;
datmean = sum(rowsums) / sum(nums_in_row);

%% replace the NaNs with values calculated from the marginals
% loop over row chunks
for i = 1:chunkdims(1):datsize(1)
    % create row chunk-selection index
    idx1 = i:min(i+chunkdims(1)-1,datsize(1));
    % loop over column chunks
    for j = 1:chunkdims(2):datsize(2)
        % construct column chunk-selection index
        idx2 = j:min(j+chunkdims(2)-1,datsize(2));
        % replace NaNs in chunk with value based on marginals
        dat = D.dat(idx1,idx2);
        [r c] = find(isnan(dat));
        x = find(isnan(dat));
        dat(x) = rowmeans(r) .* colmeans(c)' / datmean;
        D.dat(idx1,idx2) = dat;
    end
end

