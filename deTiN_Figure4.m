pairs=load_struct('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan3/pair_table.tsv') ;
samples=load_struct('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan3/sample_table.tsv');
[i m]=ismember(pairs.case_sample,samples.sample_id);
pairs.disease_type(i)=samples.primary_disease(m); pairs.disease_type=pairs.disease_type';
pairs=reorder_struct(pairs,~ismember(pairs.disease_type,{'GBM','OV'})); %less then 5 samples
boxplot(str2double(pairs.combined_deTiN_TiN_num),pairs.disease_type);

LUAD_pairs_signature_TiN=reorder_struct(pairs,(str2double(pairs.combined_deTiN_TiN_num(ismember(pairs.disease_type,'LUAD')))>.015));
LUAD_pairs_signature_TiN=reorder_struct(LUAD_pairs_signature_TiN,~ismember(LUAD_pairs_signature_TiN.maf_file_capture_deTiN,''));
M=load_struct(LUAD_pairs_signature_TiN.maf_file_capture_deTiN{1});
for i=2:slength(LUAD_pairs_signature_TiN)
m=load_struct(LUAD_pairs_signature_TiN.maf_file_capture_deTiN{i});
M=mergeStruct(M,m);
end
for i=1:slength(M)                                                     
M.change{i,1}=strcat([M.Reference_Allele{i},'->',M.Tumor_Seq_Allele2{i}]);
end
 M=reorder_struct(M,cellfun(@length,M.Tumor_Seq_Allele2)<=1);