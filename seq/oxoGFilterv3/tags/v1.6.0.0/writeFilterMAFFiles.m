function [] = writeFilterMAFFiles(mafTable, outputMAFFilename, mutationMap)
%
% mafTable must contain ALL columns that are needed for output.
% mafTable is assumed to be the annotated maf file.
% TODO: autodetect whether it is a file or struct and behave accordingly.
% TODO: Finish documentation
%

mafTable.pox = -1 * ones(length(mafTable.Start_position),1);
mafTable.qox = -1 * ones(length(mafTable.Start_position),1);
mafTable.isArtifactMode = -1 * ones(length(mafTable.Start_position),1);
mafTable.oxoGCut = -1 * ones(length(mafTable.Start_position),1);

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
    
    % Add the fields to the mafTable
    mafTable.pox(i) = newMutationInfo.pox;
    mafTable.qox(i) = newMutationInfo.qox;
    mafTable.isArtifactMode(i) = newMutationInfo.isArtifactMode;
    mafTable.oxoGCut(i) = newMutationInfo.cut;
    mafTable.i_t_Foxog(i) = newMutationInfo.foxog;
end

% TODO: Convert all (incl. struct refs above) to use a string for structure
%  fields.  That way can just provide a list of new fields and that is it.


% Write a file that contains all of the mutations 
printStruct(mafTable, 1:length(mafTable.Start_position), [outputMAFFilename '.all.maf.annotated'], mafTable.headline(1:(end-1)));

% Write a file that only contains passing mutations
passingMafTable = trimStruct(mafTable, ~mafTable.oxoGCut);
printStruct(passingMafTable, 1:length(passingMafTable.Start_position), outputMAFFilename, passingMafTable.headline(1:(end-1)));
