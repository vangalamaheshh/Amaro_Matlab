m2tsv=load_struct('/xchip/cga_home/amaro/CLL/M2_DeTiNM2/m2mafs.tsv');

m2tsv=reorder_struct(m2tsv,~ismember(m2tsv.M2_deTiN_maf,''));


indel_genes={'TP53','ATM','NOTCH1'};
for i=1:slength(m2tsv)
m2=load_struct(m2tsv.m2_maf{i});
m2D=load_struct(m2tsv.M2_deTiN_maf{i});
m2=reorder_struct(m2,ismember(m2.Variant_Type,'DEL')|ismember(m2.Variant_Type,'INS'));
m2D=reorder_struct(m2D,ismember(m2D.Variant_Type,'DEL')|ismember(m2D.Variant_Type,'INS'));

m2tsv.m2nToT(i,1)=slength(m2);
m2tsv.m2_D_nToT(i,1)=slength(m2D);
m2tsv.m2n_sig(i,1)=sum(ismember(m2.Hugo_Symbol,indel_genes));
m2tsv.m2n_D_sig(i,1)=sum(ismember(m2D.Hugo_Symbol,indel_genes));
if mod(i,5)==0
i
end
end