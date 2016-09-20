TSV=load_struct('SegFiles.tsv')
for i=1:slength(TSV)
seg=load_table(TSV.seg_path{i});
TSV.var_CNV(i,1)=var(seg.Modal_HSCN_1+seg.Modal_HSCN_2);
end
