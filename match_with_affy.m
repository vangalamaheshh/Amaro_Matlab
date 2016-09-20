function sets=match_with_affy(sets,lst)

if iscell(lst)
  lst=strvcat(lst);
end
 
for i=1:length(sets)
  sets{i}.idx=findstrings_list(lst,strvcat(sets{i}.genes));
end
