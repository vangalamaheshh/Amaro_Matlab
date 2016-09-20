
function combineValandABSmaf(abs_maf,val_maf,out_file)
abs_maf=load_struct(abs_maf);
val_maf=load_struct(val_maf);
abs_maf.key=strcat(abs_maf.Chromosome,abs_maf.Start_position);
val_maf.key=strcat(val_maf.Chromosome,val_maf.Start_position);


for i=1:slength(val_maf)
k=find(ismember(abs_maf.key,val_maf.key{i}));
if isempty(k)
val_maf.cancer_cell_fraction{i,1}='NaN';               
val_maf.ccf_CI95_low{i,1}='NaN';                                                                                                                                                     
val_maf.ccf_CI95_high{i,1}='NaN';                     
else
val_maf.cancer_cell_fraction{i,1}=abs_maf.ccf_hat{k};
val_maf.ccf_CI95_high{i,1}=abs_maf.ccf_CI95_high{k};
val_maf.ccf_CI95_low{i,1}=abs_maf.ccf_CI95_low{k};  
end
end

save_struct(val_maf,out_file) 

end

function test
%processing mafs to add ccf
abs_maf='~/Projects/CLL_Rush/ABSOLUTE/ABSOLUTE_results/CLL_Rush_4.14.14/reviewed/SEG_MAF/CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2_ABS_MAF.txt';
val_maf='/xchip/cga/gdac-prod/cga/jobResults/MutationValidator/CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2/7935165/CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2.validated.maf';
out_file='~/Projects/CLL_Rush/ABSOLUTE/CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2_combinedmaf.tsv';
combineValandABSmaf(abs_maf,val_maf,out_file)
end