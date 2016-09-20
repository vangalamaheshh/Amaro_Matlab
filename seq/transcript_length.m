function n=transcript_length(gene)
% (Gaddy's code from March2008)
res=query_spliceminer(gene);
if isempty(res)
  n=NaN;
else
  nm_idx=grep('^NM',res(:,5),1);
  n=str2num(res{nm_idx(end),10});
end
