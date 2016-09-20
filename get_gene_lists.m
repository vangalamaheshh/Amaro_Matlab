function [lsts,idx]=get_gene_lists(sets)

idx=ones(length(sets),1);
for i=1:length(sets)
  tmp=cell2mat(sets{i}.idx(~cellfun('isempty',sets{i}.idx)));
  if isempty(tmp)
    idx(i)=0;
  else
    lsts{i}=tmp;
  end
end

idx=find(idx);
lsts=lsts(idx);
