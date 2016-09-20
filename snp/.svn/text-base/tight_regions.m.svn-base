function regs=tight_regions(regs,CL21,ts)

for k=1:size(regs,1)
  for dodel=0:1
    for i=1:length(regs{k,dodel+1})
      t1=ts(k);
      if dodel
        t1=-ts(k); % log2(4/3)-1
        x=CL21.smooth((regs{k,dodel+1}(i).st):(regs{k,dodel+ ...
                            1}(i).en),:);
        x1=double(x<t1);
        pk=regs{k,dodel+1}(i).peak-regs{k,dodel+1}(i).st+1;
        keyboard
      else
        t1=ts(k);
        x=CL21.smooth((regs{k,dodel+1}(i).st):(regs{k,dodel+ ...
                            1}(i).en),:);
        x1=double(x<t1);
        pk=regs{k,dodel+1}(i).peak-regs{k,dodel+1}(i).st+1;
        keyboard
      end
    end
  end
end
