function [pairs] = retrieveUniqueCaseControlBarcodes(mafTable)
% Collect all of the tumor sample barcodes in the maf file.
%
% [pairs] = retrieveUniqueCaseControlBarcodes(mafTable)
%
% Output:
%   pairs -- structs with fields case and control, which are sample names.
%       Each field in the struct is a cell with the name.
%
%  IMPORTANT:  Assumes that the case sample is in the Tumor_Sample_Barcode
%   and the control sample is in the Matched_Norm_Sample_Barcode column.
%   Assumes that the maf file has the case/control pair on each line.
%
if length(mafTable.Tumor_Sample_Barcode) ~= length(mafTable.Matched_Norm_Sample_Barcode)
   error('Unequal number of case/control samples.') 
end

% Preallocate for speed
tmp = struct('case', cell(1), 'control', cell(1));
cases = unique(strcat(mafTable.Tumor_Sample_Barcode, ',', mafTable.Matched_Norm_Sample_Barcode ));
numCases = length(unique(strcat(mafTable.Tumor_Sample_Barcode,mafTable.Matched_Norm_Sample_Barcode)));
pairs = repmat(tmp,numCases,1);

for i = 1:numCases
    caseString = cases{i};
    tmpSplit = regexp(caseString,'(,)','split');
    pairs(i).case = tmpSplit{1};
    pairs(i).control = tmpSplit{2};
end
