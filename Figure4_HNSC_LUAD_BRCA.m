HNSC=load_struct('~/Documents/HNSC_LUAD');
BRCA=load_struct('~/Documents/BRCA');

BRCA.tumor_in_normal_estimate=str2double(BRCA.tumor_in_normal_estimate);
BRCA.TiN_LoH=str2double(BRCA.TiN_LoH);

HNSC.tumor_in_normal_estimate=str2double(HNSC.tumor_in_normal_estimate);
HNSC.combined_deTiN_TiN_num=str2double(HNSC.combined_deTiN_TiN_num);

HNSC=reorder_struct(HNSC,~(HNSC.combined_deTiN_TiN_num==-1&HNSC.tumor_in_normal_estimate==-1));
BRCA=reorder_struct(BRCA,~(BRCA.tumor_in_normal_estimate==-1));

hnsc=reorder_struct(HNSC,~cellfun(@isempty,strfind(HNSC.pair_id,'HNSC')));
luad=reorder_struct(HNSC,~cellfun(@isempty,strfind(HNSC.pair_id,'LUAD')));


pie([sum(BRCA.tumor_in_normal_estimate<.02&BRCA.tumor_in_normal_estimate>-1);sum(BRCA.tumor_in_normal_estimate>=.02&BRCA.tumor_in_normal_estimate>-1)])
pie([sum(luad.combined_deTiN_TiN_num<.02&luad.tumor_in_normal_estimate>-1);sum(luad.combined_deTiN_TiN_num>=.02)])
pie([sum(hnsc.combined_deTiN_TiN_num<.02&hnsc.tumor_in_normal_estimate>-1);sum(hnsc.combined_deTiN_TiN_num>=.02)])

boxplot([hnsc.combined_deTiN_TiN_num;luad.combined_deTiN_TiN_num;BRCA.tumor_in_normal_estimate],[ones(slength(hnsc),1);2*ones(slength(luad),1);3*ones(slength(BRCA),1)],'Labels',{'HNSC','LUAD','BRCA'})