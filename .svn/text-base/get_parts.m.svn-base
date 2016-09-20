function parts=get_parts(v,nparts)

nparts=min(nparts,length(v));
part_sz=floor(length(v)/nparts);
m=mod(length(v),nparts);
c=1;
if m~=0
  for i=1:m
    parts{i}=v(c:(c+part_sz));
    c=c+part_sz+1;
  end
  for i=(m+1):nparts
    parts{i}=v(c:(c+part_sz-1));
    c=c+part_sz;    
  end
else
  for i=1:nparts
    parts{i}=v((1:part_sz)+(i-1)*(part_sz));
  end
end

