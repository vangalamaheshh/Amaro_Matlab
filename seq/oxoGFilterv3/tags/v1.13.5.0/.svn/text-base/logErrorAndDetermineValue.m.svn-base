function newValue = logErrorAndDetermineValue(mafTable, i, mutationInfo, mafTableFieldName, mutationInfoFieldName, errorFid, NULL_VALUE)
% Returns the value to use in the mafTable
% errorFid whould already be open in write mode.
%
% TODO: Finish documentation

mafTableVals = getfield(mafTable, mafTableFieldName);
originalValue = mafTableVals(i);
newValue = getfield(mutationInfo, mutationInfoFieldName);

if (originalValue ~= NULL_VALUE) && (~strcmp(num2str(originalValue), num2str(newValue)))
   if iscell(mafTable.Chromosome)
       chr =  mafTable.Chromosome{i};
   else
       chr =  mafTable.Chromosome(i);
   end
   fprintf(errorFid, 'WARNING: Prexisting value being changed in the MAF file for %s:%d (line %d): %s %f --> %f\n', chr, mafTable.Start_position(i), i, mafTableFieldName, originalValue,  newValue );
end
