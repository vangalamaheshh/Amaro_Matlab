function [] = writeCaseMafTempFile(mafTable, unameIndices, M, fdrStruct, acVal, acs, outputCaseTempFilename)
%
% TODO: Finish Documentation
%
tmpMafTable = trimStruct(mafTable, unameIndices);
tmpMafTable.pox = fdrStruct.pox;
tmpMafTable.qox = fdrStruct.qox;
tmpMafTable.cut = fdrStruct.cut;
tmpMafTable.isArtifactMode = M.isArtifactMode;
    % TODO: Still need pox_cutoff
    
printStruct(tmpMafTable, 1:length(tmpMafTable.pox), outputCaseTempFilename);
