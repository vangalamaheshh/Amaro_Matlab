function [x,y]=remove_large_segments(C,x,sz)

y=x;
for j=1:max(C.chrn)
  inchr=find(C.chrn==j);
  if ~isempty(inchr)
    for i=1:size(x,2)
      rl=runlength(x(inchr,i));
      rlsz=rl(:,2)-rl(:,1)+1;
      large=find(rlsz>=sz & rl(:,3)>0);
      if ~isempty(large)
        rl(large,3)=0;
      end
      x(inchr,i)=derunlength(rl);
    end
    disp([ i j]);
  end
end
y=y-x;
