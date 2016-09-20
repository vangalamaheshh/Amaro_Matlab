M_ABS=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/ABS_Mafs/AggregatedGCT.maf');
MSIG_Pairs=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/MSIGfirehosetable.tsv');
ClinicalData=load_struct('GCT_Armlvl_data/ClinicalPathData.txt');
All_Pairs=load_struct('~/Documents/All_Pairs_All_Figures.txt');
M_ABS=reorder_struct(M_ABS,ismember(M_ABS.sample_name,All_Pairs.case_sample));
for i=1:slength(M_ABS)
M_ABS.bsp_id{i,1}=M_ABS.sample_name{i}(end-7:end);
end
M_ABS=reorder_struct(M_ABS,ismember(M_ABS.i_failure_reasons,{'NA',''}));
M_ABS.alt=str2double(M_ABS.alt);
M_ABS.ref=str2double(M_ABS.ref);
M_ABS.tumor_f=M_ABS.alt./(M_ABS.alt+M_ABS.ref);
M_ABS=reorder_struct(M_ABS,M_ABS.tumor_f>.05);
M_ABS_msig=reorder_struct(M_ABS,ismember(M_ABS.sample_name,MSIG_Pairs.case_sample)); % for mutsig plot
location=load_struct('LocationFigureForSamples.txt');
[table.count,table.id]=count(M_ABS_msig.bsp_id);

for i=1:slength(table)
if ~isempty(find(ismember(table.id{i},location.uniq_samp_tcga_barcode)));
table.location{i,1}=location.location_testes{ismember(location.uniq_samp_tcga_barcode,table.id{i})};
end
end
table.location(38:49)={'yes'};
boxplot(table.count./28.7,table.location)
title('Location of Sample Sequenced','FontSize',20)
set(gca,'XTickLabel',{'Other','Testes'},'FontSize',20)

ABS_Table=load_struct('/Volumes/xchip_cga_home/amaro/TestesRotation/FreezeSet/AggregatedGCTABSTable.txt');

res=regexp(ABS_Table.array,'SM-[0-9][A-Z]+[0-9]*','match');

for i=1:length(res)
if length(res{i})>1
ABS_Table.sm_id{i,1}=res{i}{1};
else 
ABS_Table.sm_id{i,1}='NA';
end
end

for i=1:slength(table)
table.purity(i,1)=str2double(ABS_Table.purity(ismember(ABS_Table.sm_id,table.id{i})));
end
for i=1:slength(table)
k=find(ismember(ABS_Table.sm_id,table.id{i}));
table.pair_id{i,1}=ABS_Table.array{k};
end
for i=1:slength(table)
    table.ind_id(i,1)=regexp(table.pair_id{i},'DFCI_[0-9]+','match');
end
for i=1:slength(table)
    k=find(ismember(ClinicalData.FH_names,strcat(['Testes_',table.ind_id{i}])));
    table.vital_status(i,1)=ClinicalData.Vitalstatus(k);
end






table.score(ismember(table.location,'yes')&ismember(table.vital_status,'1'))=1;
table.score(ismember(table.location,'yes')&ismember(table.vital_status,'2'))=2;
table.score(ismember(table.location,'no')&ismember(table.vital_status,'1'))=3;
table.score(ismember(table.location,'no')&ismember(table.vital_status,'2'))=4;


ABS_SEG_primaries=reorder_struct(ABS_SEG,ismember(ABS_SEG.pair_id,table.pair_id(ismember(table.location,'yes'))));
ABS_SEG_non_testes=reorder_struct(ABS_SEG,ismember(ABS_SEG.pair_id,table.pair_id(~ismember(table.location,'yes'))));

P53_loc=[7569720,7592868];
P53_loc(3)=P53_loc(2)-(P53_loc(2)-P53_loc(1));
MDM2_loc=[69199952,69241324];
MDM2_loc(3)=MDM2_loc(2)-(MDM2_loc(2)-MDM2_loc(1));

primaries=unique(ABS_SEG_primaries.pair_id);
non_testes=unique(ABS_SEG_non_testes.pair_id);
for i=1:length(primaries)
    s=reorder_struct(ABS_SEG_primaries,ismember(ABS_SEG_primaries.pair_id,primaries{i}));
    CN_tp53(i,1)=s.rescaled_total_cn((s.Chromosome==17&s.Startbp<P53_loc(1)&s.Endbp>P53_loc(2))|(s.Chromosome==17&s.Startbp<P53_loc(3)&s.Endbp>P53_loc(3)));
    if isempty(s.rescaled_total_cn((s.Chromosome==12&s.Startbp<MDM2_loc(1)&s.Endbp>MDM2_loc(2))|(s.Chromosome==12&s.Startbp<MDM2_loc(3)&s.Endbp>MDM2_loc(3))))
        CN_MDM2(i,1)=NaN;
    else
    CN_MDM2(i,1)=s.rescaled_total_cn((s.Chromosome==12&s.Startbp<MDM2_loc(1)&s.Endbp>MDM2_loc(2))|(s.Chromosome==12&s.Startbp<MDM2_loc(3)&s.Endbp>MDM2_loc(3)));
    end
end





for i=1:length(non_testes)
    s=reorder_struct(ABS_SEG_non_testes,ismember(ABS_SEG_non_testes.pair_id,non_testes{i}));
    if isempty(s.rescaled_total_cn((s.Chromosome==17&s.Startbp<P53_loc(1)&s.Endbp>P53_loc(2))|(s.Chromosome==17&s.Startbp<P53_loc(3)&s.Endbp>P53_loc(3))))
        CN_tp53(i,1)=2;
    else
    CN_tp53(i,1)=s.rescaled_total_cn((s.Chromosome==17&s.Startbp<P53_loc(1)&s.Endbp>P53_loc(2))|(s.Chromosome==17&s.Startbp<P53_loc(3)&s.Endbp>P53_loc(3)));
    end
    if isempty(s.rescaled_total_cn((s.Chromosome==12&s.Startbp<MDM2_loc(1)&s.Endbp>MDM2_loc(2))|(s.Chromosome==12&s.Startbp<MDM2_loc(3)&s.Endbp>MDM2_loc(3))))
        CN_MDM2(i,1)=2;
    else
    CN_MDM2(i,1)=s.rescaled_total_cn((s.Chromosome==12&s.Startbp<MDM2_loc(1)&s.Endbp>MDM2_loc(2))|(s.Chromosome==12&s.Startbp<MDM2_loc(3)&s.Endbp>MDM2_loc(3)));
    end
end




