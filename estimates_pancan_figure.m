pancan = load_table('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan_NT_paper/TiN_Estimates/estimates.tsv');
breast = load_table('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan_NT_paper/TiN_Estimates/estimates.brca.tsv');
pancan=rmfield(pancan,{'header','headline'}); breast=rmfield(breast,{'header','headline'}); 

pancan=reorder_struct(pancan,~(pancan.tumor_in_normal_estimate==-1&ismember(pancan.normal_tissue,'NB')));

for i=1:slength(pancan)
pancan.normal_tissue{i,1}=pancan.pair_id{i}(end-1:end);
strs=regexp(pancan.pair_id{i},'_','split');
pancan.tumor_type{i,1}=strs{3};
pancan.tumor_sample{i,1}=strs{1};
end

breast.normal_tissue=repmat({'NT'},slength(breast),1);
breast.tumor_type=repmat({'BRCA'},slength(breast),1);
figure()
pie([sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NT'))<.02)+sum(breast.combined_deTiN_TiN_num<.02),sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NT'))>=.02)+sum(breast.combined_deTiN_TiN_num>=.02)])
figure()
pie([sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NB'))<.02),sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NB'))>=.02)])


tumor_types=unique(pancan.tumor_type);
for i = 1:length(tumor_types)
    figure()
    pie([sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NT')&ismember(pancan.tumor_type,tumor_types{i}))<.02),sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NT')&ismember(pancan.tumor_type,tumor_types{i}))>=.02)],[0,0],num2cellstr([sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NT')&ismember(pancan.tumor_type,tumor_types{i}))<.02),sum(pancan.combined_deTiN_TiN_num(ismember(pancan.normal_tissue,'NT')&ismember(pancan.tumor_type,tumor_types{i}))>=.02)]))
    title(tumor_types{i})
end

pie([sum(breast.combined_deTiN_TiN_num<.02),sum(breast.combined_deTiN_TiN_num>=.02)])


[n,tsamples]=count(pancan.tumor_sample);
trio_tumors=tsamples(n>1);
find(ismember(pancan.tumor_sample,trio_tumors)&pancan.combined_deTiN_TiN_num>=.02)


