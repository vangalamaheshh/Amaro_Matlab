function d=segment_dist(C,x)

d=cell(size(x,2),1);
for j=1:max(C.chrn)
  inchr=find(C.chrn==j);
  if ~isempty(inchr)
    for i=1:size(x,2)
      rl=runlength(x(inchr,i));
      rlsz=rl(:,2)-rl(:,1)+1;
      d{i}=[ d{i}; rlsz(find(rl(:,3)>0)) ];
    end
  end
end
