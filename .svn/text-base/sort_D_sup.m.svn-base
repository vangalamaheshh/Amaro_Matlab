function cidx=sort_D_sup(D,supid,types,by_name)

if length(supid)==1
  types={types};
end

if exist('by_name','var')
  if iscell(D.sdesc)
    [tmp,tidx]=sortrows(strvcat(D.sdesc));
  else
    [tmp,tidx]=sortrows(D.sdesc);
  end
  D=reorder_D_cols(D,tidx);
else
  tidx=1:size(D.dat,2);
end

cidx=1:size(D.dat,2);

for j=length(supid):-1:1
  t=types{j};
  s=supid(j);
  if isempty(t)
    t=unique_keepord(D.supdat(s,:));
  end
  idx=[];
  for i=1:length(t)
    idx=[idx find(D.supdat(s,cidx)==t(i))];
  end
  cidx=cidx(idx);
end

cidx=tidx(cidx);

