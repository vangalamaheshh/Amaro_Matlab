NewABS=load_struct('/Users/amaro/Documents/CLL_Questions_For_Dan/Stilgenbauerfullset.1.17.ATW.ABSOLUTE.Table.txt');
OldABS=load_struct('/Users/amaro/Documents/CLL_Questions_For_Dan/Stilgenbauerfullset.11.21.ATW.ABSOLUTE.tableTP2Later.txt');
NewABS=load_table('/Users/amaro/Documents/CLL_Questions_For_Dan/Stilgenbauerfullset.1.17.ATW.ABSOLUTE.Table.txt');
OldABS=load_table('/Users/amaro/Documents/CLL_Questions_For_Dan/Stilgenbauerfullset.11.21.ATW.ABSOLUTE.tableTP2Later.txt');
slength(NewABS)
slength(OldABS)
NewABS=reorder_struct(ismember(NewABS.array,OldABS.array));
NewABS=reorder_struct(NewABS.array,ismember(NewABS.array,OldABS.array));
NewABS=reorder_struct(NewABS,ismember(NewABS.array,OldABS.array));
slength(NewABS)
OldABS=reorder_struct(OldABS,ismember(OldABS.array,NewABS.array));
slength(OldABS)
[i n]=ismember(NewABS.array,OldABS.array);
i
n
figure()
plot(OldABS.purity,NewABS.purity,'b.')
xlabel('Old Purity','FontSize',18)
ylabel('New Purity','FontSize',18)
NewABS.Purity_diff=abs(OldABS.purity-NewABS.purity);
hist(NewABS.Purity_diff)
hist(NewABS.Purity_diff,30)
whitelist_maf=load_struct('/Users/amaro/Documents/CLL_Questions_For_Dan/An_CLL_April_v1.final_analysis_set.maf.txt');
whitelist_maf
whitelist_maf=reorder_struct(whitelist_maf,ismember(whitelist_maf.Hugo_Symbol,{'SF3B1','NOTCH1','XPO1','RPS15','ZMYM3','CHD2'}))
whitelist_maf
count(whitelist_maf.Tumor_Sample_Barcode)
whitelist_maf_sub=reorder_struct(whitelist_maf,ismember(whitelist_maf.Tumor_Sample_Barcode,{'CW145-Tumor','CLL_166-Tumor','CW147-Tumor','CW230-Tumor','CW238-Tumor','CLL_170-Tumor','CLL_101-Tumor','CLL_103-Tumor'}));

for i=24:length(GermlineMafs.pair_id)
if exist(GermlineMafs.germline_maf_singlesample_analysisready{i},'file')
g=load_struct(GermlineMafs.germline_maf_singlesample_analysisready{i});
GermlineMafs.HetsOnX(i,1)=sum(ismember(g.genotype,'0/1')&ismember(g.Chromosome,'X'));
GermlineMafs.HomsOnX(i,1)=sum(ismember(g.genotype,'1/1')&ismember(g.Chromosome,'X'));
i
end
end