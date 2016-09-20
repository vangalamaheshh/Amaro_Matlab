function [idx,had_a_match]=get_D_idx(D,dim,lst)

had_a_match=[];
idx=[];
switch dim
 case {'col','cols','samples','conditions'}
  strs=D.sdesc;
 case {'row','rows','genes','features'}
  strs=D.gacc;
end

if iscell(lst)
  nlst=length(lst);
else
  nlst=size(lst,1);
end

for i=1:nlst
  if iscell(lst)
    lsti=deblank(lst{i});
  else
    lsti=deblank(lst(i,:));
  end
  pos=strmatch(lsti,strs,'exact');
  if ~isempty(pos)
    if length(pos)>1
      assert(['found ' lsti ' more than once']);
    else
      idx=[idx pos];
      had_a_match=[had_a_match i];
    end
  end
end


