function write_regs_calls(CL,fsuffix,regs,calls,pvs)


for k=1:2
  les_acc=[];
  for i=1:length(regs{k})
    st=regs{k}(i).st;
    en=regs{k}(i).en;
    pk=regs{k}(i).peak;
    pkst=sprintf('%g-%g ',CL.pos(regs{k}(i).peak_st),CL.pos(regs{k}(i).peak_en));
    les_st=sprintf('%s:%g Mb - %s:%g Mb  peak at %s:%s, fdr %f : ( %d:%d )',...
                   CL.chr{st},CL.pos(st),CL.chr{en},CL.pos(en),CL.chr{pk}, ...
                   pkst,pvs{k}(regs{k}(i).peak),regs{k}(i).peak_st,regs{k}(i).peak_en); 
    les_acc{i}=les_st;
  end
  if k>1
    write_eisen_dat(['deletions_log_' num2str(size(calls{k},2)) '_' fsuffix '.txt'],...
                    strvcat(les_acc),repmat(' ',length(regs{k}),1),strvcat(CL.sdesc),'Deletions', ...
                    calls{k},[],[],1);
  else
    write_eisen_dat(['amplifications_log_' num2str(size(calls{k},2)) '_' fsuffix '.txt'],...
                    strvcat(les_acc),repmat(' ',length(regs{k}),1),strvcat(CL.sdesc),'Amplifications', ...
                    calls{k},[],[],1);
  end
end
