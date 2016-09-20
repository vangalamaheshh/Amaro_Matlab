function Fi = pick_highest_priority_mutation(F,type_priority)
% pick_highest_priority_mutation(F,type_priority)
% default type_priority = {'Nonsense','Read-through','Missense','Splice_site',...
%    'Synonymous','5''-UTR','3''-UTR','Intron'}

if ~exist('type_priority','var')
  type_priority = {'Nonsense','Read-through','Missense','Splice_site',...
    'Synonymous','5''-UTR','3''-UTR','Intron'};
end

for p=1:length(type_priority)
  idx = find(strcmpi(F.type,type_priority{p}),1);
  if ~isempty(idx)
    Fi = reorder_struct(F,idx);
    break;
end,end

