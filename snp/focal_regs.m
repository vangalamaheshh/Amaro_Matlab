function calls=focal_regs(regs,)

for k=1:size(regs,1)
  for dodel=0:(size(regs,2)-1)
    for i=1:length(regs{k,dodel+1})
      t1=ts(k);
      if dodel
        t1=-ts(k); % log2(4/3)-1
        calls{k,dodel+1}(i,:)=double(CL21.smooth(regs{k,dodel+1}(i).peak,:)<t1);
      else
        t1=ts(k);
%        keyboard
        calls{k,dodel+1}(i,:)=double(CL21.smooth(regs{k,dodel+1}(i).peak,:)>t1);
      end
    end
  end
end
