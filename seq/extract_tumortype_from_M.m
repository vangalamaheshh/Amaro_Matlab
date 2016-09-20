function G=extract_tumortype_from_M(M,tti)
  G=[]; G.rank=(1:slength(M.gene))'; G.name = M.gene.name;
  G.nnon = sum(M.ttype.gene_nmut_non(tti,:,end),1)'; G.Nnon = round(sum(M.ttype.gene_non_repcov(tti,:,end),1))';
  G.nsil = sum(M.ttype.gene_nmut_sil(tti,:,end),1)'; G.Nsil = round(sum(M.ttype.gene_sil_repcov(tti,:,end),1))';
  G.nflank = sum(M.ttype.gene_nmut_flank(tti,:,end),1)'; G.Nflank = round(sum(M.ttype.gene_flank_repcov(tti,:,end),1))';
