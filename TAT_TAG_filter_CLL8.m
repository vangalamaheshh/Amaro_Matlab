aM=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/LEGO_Plot_Clonal/MutsigOut/final_analysis_set.maf');
%M.key=strcat(M.Start_Position,M.Tumor_Sample_Barcode);
affected_tumors={'CLL-CW50-Tumor-SM-UFPS','CLL-CW53-Tumor-SM-UFPY','CLL-CW54-Tumor-SM-UFQ1',...
    'CLL-CW56-Tumor-SM-UFQ5','CLL-CW62-Tumor-SM-UFRK','CLL-CW73-Tumor-SM-UFS7','CLL-CW95-Tumor-SM-UFTG',...
    'CLL-CW103-Tumor-SM-UFTW','CLL-CW107-Tumor-SM-Z133','CLL-CW109-Tumor-SM-Z137','CLL-CW80-Tumor-SM-UFSL','CLL-CW54-Tumor-SM-UFQ1',...
'CLL-CW48-Tumor-SM-UFPO'};
%M=reorder_struct(M,ismember(M.Variant_Type,'SNP'));
%M.mut_base=strcat(M.Reference_Allele,M.Tumor_Seq_Allele2);
contexts={'15','16','53','49'};
alt={'A','T'};
aM=reorder_struct(aM,(ismember(aM.Tumor_Seq_Allele2,alt)&ismember(aM.context65,contexts)&ismember(aM.patient,affected_tumors)&ismember(aM.Variant_Type,'SNP')));


aM.key=strcat(aM.Tumor_Sample_Barcode,aM.Start_Position);

M=reorder_struct(M,~ismember(M.key,aM.key));


for i=1:slength(aM)
    if isequal(aM.Reference_Allele{i},'A')
        aM.next_base{i,1}=aM.ref_context{i}(12);
    else
        aM.next_base{i,1}=aM.ref_context{i}(10);
    end
end