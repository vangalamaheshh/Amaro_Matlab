AddAnnot=load_struct('/Users/amaro/Downloads/MAF_fileJMML.tsv');
J_maf=load_struct('/Users/amaro/Dropbox/Broad Shared/All_Pairs.final_analysis_set.filtered.maf');
J_maf.key=strcat(J_maf.Tumor_Sample_Barcode,J_maf.Start_position);
AddAnnot.key=strcat(AddAnnot.Tumor_Sample_Barcode,AddAnnot.Start_position);
