function C=join_close(C,min_delta,t)

for c=1:23
  ci{c}=find(C.chrn==c);
end
chrs=as_row(sort(unique(C.chrn)));
for i=1:size(C.cbs,2)
  for c=chrs
    rl=runlength(C.cbs(ci{c},i)');
    if size(rl,1)>1
      update=0;
      
      delta=diff(rl(:,3));
      above_t=(rl(:,3)>t);
      both_above_t=above_t(2:end) & above_t(1:(end-1));
      join_seg=both_above_t & (abs(delta)<=min_delta);
      if any(join_seg)
        rl_join=runlength(join_seg');
        for k=find(rl_join(:,3)==1)'
          [i c k 1]
          rli=rl_join(k,1):(rl_join(k,2)+1);
          sz=rl(rli,2)-rl(rli,1)+1;
          rl(rli,3)=(rl(rli,3)'*sz)/sum(sz);
        end
        update=1;
      end
      
      above_t=(rl(:,3)<-t);
      both_above_t=above_t(2:end) & above_t(1:(end-1));
      join_seg=both_above_t & (abs(delta)<=min_delta);
      if any(join_seg)
        rl_join=runlength(join_seg');
        for k=find(rl_join(:,3)==1)'
          [i c k -1]
          rli=rl_join(k,1):(rl_join(k,2)+1);
          sz=rl(rli,2)-rl(rli,1)+1;
          rl(rli,3)=(rl(rli,3)'*sz)/sum(sz);
        end
        update=1;
      end
      
      if update
        C.cbs(ci{c},i)=(derunlength(rl))';
        disp('x');
      end
    end
  end
end
