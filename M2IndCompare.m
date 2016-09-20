M=load_struct('~/Projects/CLL/Indel_M2_Compare/IndelLM2.tsv')

for i=1:slength(M)

m2=load_struct(M.m2_maf{i});
m2=reorder_struct(m2,ismember(m2.Variant_Classification,'Frame_Shift_Ins')|ismember(m2.Variant_Classification,'Frame_Shift_Del')|ismember(m2.Variant_Classification,'In_Frame_Del')|ismember(m2.Variant_Classification,'In_Frame_Ins'));
ind=load_struct(M.indel_maf_file_capture{i});
ind=reorder_struct(ind,ismember(ind.Variant_Classification,'Frame_Shift_Ins')|ismember(ind.Variant_Classification,'Frame_Shift_Del')|ismember(ind.Variant_Classification,'In_Frame_Del')|ismember(ind.Variant_Classification,'In_Frame_Ins'));
M.m2_ATM_TP53_NOTCH1_count(i,1)=sum(ismember(m2.Hugo_Symbol,sig_genes));
M.ind_ATM_TP53_NOTCH1_count(i,1)=sum(ismember(ind.Hugo_Symbol,sig_genes));
M.m2Total(i,1)=slength(m2);
M.indTotal(i,1)=slength(ind);
if mod(i,5)==0
    i
end
end