function [mafTable] = writeFilterMAFFiles(mafTable, outputMAFFilename, mutationMap)
%
% [mafTable] = writeFilterMAFFiles(mafTable, outputMAFFilename, mutationMap)
%   Responsible for two actions:
%       1) write filtered and unfiltered maf files given a maf struct
%       created by loadMAFTable.
%       2) Returns a new mafTable that includes all of the original
%       mafTable, plus new fields (pox, qox, isArtifactMode, oxoGCut,
%       i_t_Foxog, pox_cutoff).
%
% mafTable must contain ALL columns that are needed for output.
% 
% TODO: autodetect whether it is a file or struct and behave accordingly.
% 
% Writes any warnings or errors to [outputMAFFilename].err.txt
%
% See also load_table, loadMAFTable

errorFid = fopen([outputMAFFilename '.err.txt' ], 'w');

% Null value should be a value that would never appear in the mutation map
%   fields
NULL_VALUE = -1;

% List of fields that are going to be added to the output maf.
fieldList = [{'pox'};{'qox'};{'pox_cutoff'};{'isArtifactMode'};{'oxoGCut'};];
for i = 1:length(fieldList)
   fieldName = fieldList{i};
   
   % Field does not already exist, so initialize it to a value that would 
   %    never already be in the mafTable.
   if ~isfield(mafTable, fieldName)
       mafTable = setfield(mafTable, fieldName, NULL_VALUE * ones(length(mafTable.Start_position),1));
   end
end


% Go through every mutation in the given maf file.
for i = 1:length(mafTable.Start_position)
    
    % Acquire the key into the mutation map for this mutation
    if iscell(mafTable.Chromosome)
        k = createMutationMapKey(mafTable.Tumor_Sample_Barcode{i}, mafTable.Matched_Norm_Sample_Barcode{i}, mafTable.Chromosome{i}, num2str(mafTable.Start_position(i)));
    else
        k = createMutationMapKey(mafTable.Tumor_Sample_Barcode{i}, mafTable.Matched_Norm_Sample_Barcode{i}, mafTable.Chromosome(i), num2str(mafTable.Start_position(i)));
    end
    
    % Get the additional mutation information for this case and mutation.
    newMutationInfo = mutationMap(k);
    
    % Check to see if mafTable already has a value.  If it differs from the
    %   new value, write a warning.    
    % Add the fields to the mafTable    
    mafTable.pox(i) = logErrorAndDetermineValue(mafTable, i, newMutationInfo, 'pox', 'pox', errorFid, NULL_VALUE);
    mafTable.qox(i) = logErrorAndDetermineValue(mafTable, i, newMutationInfo, 'qox', 'qox', errorFid, NULL_VALUE);
    mafTable.isArtifactMode(i) = logErrorAndDetermineValue(mafTable, i, newMutationInfo, 'isArtifactMode', 'isArtifactMode', errorFid, NULL_VALUE);
    mafTable.oxoGCut(i) = logErrorAndDetermineValue(mafTable, i, newMutationInfo, 'oxoGCut', 'cut', errorFid, NULL_VALUE);
    mafTable.i_t_Foxog(i) =  logErrorAndDetermineValue(mafTable, i, newMutationInfo, 'i_t_Foxog', 'foxog', errorFid, NULL_VALUE);
    mafTable.pox_cutoff(i) = logErrorAndDetermineValue(mafTable, i, newMutationInfo, 'pox_cutoff', 'pox_cutoff', errorFid, NULL_VALUE);
end

% Write a file that contains all of the mutations 
mafTable.headline{end+1}=mafTable.headline{end};
mafTable.headline{end-1}='## OxoG Filter v3';
printStruct(mafTable, 1:length(mafTable.Start_position), [outputMAFFilename '.all.maf.annotated'], mafTable.headline(1:(end-1)));

% Write a file that only contains passing mutations
passingMafTable = trimStruct(mafTable, ~mafTable.oxoGCut);

passingMafTable.headline=mafTable.headline;

passingMafTable.headline{end+1}=passingMafTable.headline{end};
passingMafTable.headline{end-1}='## OxoG Filter v3';

printStruct(passingMafTable, 1:length(passingMafTable.Start_position), outputMAFFilename, passingMafTable.headline(1:(end-1)));

fclose(errorFid);
