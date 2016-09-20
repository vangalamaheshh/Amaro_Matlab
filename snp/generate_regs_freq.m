function regs=generate_regs_freq(CL,q,qv_thresh)
regs={[],[]};
for k=1:2
  for i=1:max(CL.chrn)
    in_chr=find(CL.chrn==i);
    qc=q{k}(in_chr);
    rl=runlength(qc<=qv_thresh);
    
    for j=(find(rl(:,3)==1)')
      mn=min(qc(rl(j,1):rl(j,2)));
      rlm=runlength(qc(rl(j,1):rl(j,2))==mn);
      for n=(find(rlm(:,3)==1)')
        [k i j n]
        regs{k}(end+1).st=in_chr(1)-1+rl(j,1);
        regs{k}(end).en=in_chr(1)-1+rl(j,2);
        regs{k}(end).peak_st=in_chr(1)-1+rl(j,1)-1+rlm(n,1);
        regs{k}(end).peak_en=in_chr(1)-1+rl(j,1)-1+rlm(n,2);
      end
    end
  end
end
