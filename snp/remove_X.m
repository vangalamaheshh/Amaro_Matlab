function C=remove_X(C)

C=reorder_D_rows(C,C.chrn~=23);
for i=1:size(C.dat,2)
  tmp=C.cbs_rl{i};
  if size(tmp,2)~=4
    error('must have chrmosome in cbs_rl');
  end
  C.cbs_rl{i}=tmp(tmp(:,4)~=23,:);
end
