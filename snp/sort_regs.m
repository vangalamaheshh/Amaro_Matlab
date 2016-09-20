function regs_sort=sort_regs(regs)

regs_sort=regs;
for k=1:2
  if ~isempty(regs{k})
    pos=cat(1,regs{k}.peak);
    [spos,sposi]=sort(pos);
    regs_sort{k}=regs_sort{k}(sposi);
  end
end
