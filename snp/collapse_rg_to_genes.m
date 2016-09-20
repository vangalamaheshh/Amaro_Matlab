function rg1=collapse_rg_to_genes(rg)

s={rg.symb};
[u,ui,uj]=unique(strvcat(s),'rows');
rg1=rg(ui);
for i=1:length(ui)
  if mod(i,1000)==0
    disp(i);
  end
  idx=find(uj==i);
  if length(idx)>1
    st=min(cat(1,rg(idx).start));
    en=min(cat(1,rg(idx).end));
    cds_st=min(cat(1,rg(idx).cds_start));
    cds_en=min(cat(1,rg(idx).cds_end));
    refseq=list2str({rg(idx).refseq},{'///'});
    status=list2str({rg(idx).status},{'///'});
    rg1(i).start=st;
    rg1(i).end=en;
    rg1(i).cds_start=cds_st;
    rg1(i).cds_end=cds_en;
    rg1(i).refseq=refseq;
    rg1(i).status=status;
  end
end

    
