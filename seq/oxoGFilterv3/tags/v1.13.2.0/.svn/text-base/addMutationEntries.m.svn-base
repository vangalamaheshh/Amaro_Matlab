function [mutMap] = addMutationEntries( pair, unameMafTable, M, fdrStruct, acVal, acs, PoxoG)
%
% Updates a structure to aggregate all data used in the maf output.
% 
%
% Key is [case '_' control '_' chromosome '_' startPosition]
% Output is a dictionary.

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
   
   % Determine the estimated pox cutline
   if all(acs < M.alt_read_count(i))
      acIndex = length(acs);
   elseif all(acs > M.alt_read_count(i))
       acIndex = 1;
   else
      acIndex = (acs == M.alt_read_count(i));
   end
   value.pox_cutoff = binocdf(acVal(acIndex), acs(acIndex), PoxoG);
   
   mutMap(keyString) = value;
end