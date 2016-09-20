function [mutMap] = addMutationEntries( pair, unameMafTable, M, fdrStruct, acVal, acs)
%
% TODO: Finish Documentation
% Updates a structure to aggregate all data used in the maf output.
% 
%
% Key is [case '_' control '_' chromosome '_' startPosition]
% Output is a dictionary.
    % TODO: Still need pox_cutoff

mutMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:length(fdrStruct.pox)

   if iscell(unameMafTable.Chromosome)
       keyString = createMutationMapKey(pair.case, pair.control, unameMafTable.Chromosome{i}, num2str(unameMafTable.Start_position(i)));
   else
       keyString = createMutationMapKey(pair.case, pair.control, unameMafTable.Chromosome(i), num2str(unameMafTable.Start_position(i)));
   end
   value.pox = fdrStruct.pox(i);
   value.qox = fdrStruct.qox(i);
   value.cut = fdrStruct.cut(i);
   value.isArtifactMode = M.isArtifactMode(i);
   value.foxog = M.foxog(i);
%    value.pox_cutoff = 
   mutMap(keyString) = value;
end