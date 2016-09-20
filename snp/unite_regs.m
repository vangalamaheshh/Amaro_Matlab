function regs=unite_regs(r)

regs=r{1};
for i=2:length(r)
  for k=1:2
    regs{k}=[regs{k} r{i}{k}];
  end
end
