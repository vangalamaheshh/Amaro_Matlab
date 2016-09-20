function [ pair_maf_list ] = split_aggregated_maf( MAF, field, list_1, list_2 )
%Splits a maf file into smaller mafs by one field
%returns a list of the maf files with their corresponding id found in the
%specified field 
% list 1 and list 2 specify the headers of the list 


MAF=load_struct(MAF);
splits=unique(MAF.(field));
mkdir('./Split_mafs');

for i=1:length(splits)
    maf=trimStruct(MAF,ismember(MAF.(field),splits{i})); maf=rmfield(maf,'N');
    save_struct(maf,sprintf('%s/Split_mafs/%s.split.maf',pwd,splits{i}));
    pair_maf_list.(list_1){i,1}=splits{i};
    pair_maf_list.(list_2){i,1}=sprintf('%s/Split_mafs/%s.split.maf',pwd,splits{i});
end


save_struct(pair_maf_list,'pair_maf_list.txt')




end



function test 

MAF=('/Users/amaro/Downloads/All_TM_NT_TP_NT.indel.strelka.maf.annotated');

pair_mafs=split_aggregated_maf(MAF,'Tumor_Sample_Barcode','pair_id','Strelka_Indel_Maf'); % Tumor Sample Barcode contains the unique identifiers used to split the maf
    
end