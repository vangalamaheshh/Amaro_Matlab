function C=recalc_rl(C)

C.cbs_rl=runlength(C.dat);
for i=1:length(C.cbs_rl)
  C.cbs_rl{i}(:,4)=C.chrn(C.cbs_rl{i}(:,1));
end
