function [regs1_matched,regs2_matched,regs_only_1,regs_only_2]=match_regs(CL21,regs1,regs2)

in_regs1=cell(1,2);
for k=1:2
  in_regs1{k}=zeros(size(CL21.dat,1),1);
  for i=1:length(regs1{k})
    in_regs1{k}(regs1{k}(i).peak_wide_st:regs1{k}(i).peak_wide_en)=i;
  end
end

in_regs2=cell(1,2);
for k=1:2
  in_regs2{k}=zeros(size(CL21.dat,1),1);
  for i=1:length(regs2{k})
    in_regs2{k}(regs2{k}(i).peak_wide_st:regs2{k}(i).peak_wide_en)=i;
  end
end



regs_only_1=cell(1,2);
regs1_matched=cell(1,2);
for k=1:2
  for i=1:length(regs1{k})
    if ~any(in_regs2{k}(regs1{k}(i).peak_wide_st:regs1{k}(i).peak_wide_en))
      regs_only_1{k}(end+1)=regs1{k}(i);
    else
      regs1_matched{k}(end+1)=regs1{k}(i);
%      disp([k i regs_H{k}(i).chrn regs_H{k}(i).peak_st unique(in_regs{k}(regs_H{k}(i).peak_wide_st:regs_H{k}(i).peak_wide_en))']);
    end
  end
end

regs_only_2=cell(1,2);
regs2_matched=cell(1,2);
for k=1:2
  for i=1:length(regs2{k})
    if ~any(in_regs1{k}(regs2{k}(i).peak_wide_st:regs2{k}(i).peak_wide_en))
      regs_only_2{k}(end+1)=regs2{k}(i);
    else
      regs2_matched{k}(end+1)=regs2{k}(i);
    end
  end
end

