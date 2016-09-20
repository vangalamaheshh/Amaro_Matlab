function res=compare_regs(r1,r2)
tmp=[];
for k=1:2
  for i=1:length(r2{k})
    tmp=[ tmp; i k r2{k}(i).peak_wide_st r2{k}(i).peak_wide_en];
  end
end

tmp1=[];
for k=1:2
  for i=1:length(r1{k})
    tmp1=[ tmp1; i k r1{k}(i).peak_wide_st r1{k}(i).peak_wide_en];
  end
end

res=tmp1-tmp;
