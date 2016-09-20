function afplots(maf,order,purity_values,purity_err)

subplot(2,1,1)
barwitherr(purity_err-purity_values,purity_values)
set(gca,'XTickLabel','')
set(gca,'YLim',[0,1]);
ylabel('Purity')

subplot(2,1,2)
boxplot(maf.i_tumor_f,maf.patient,'labelorientation','inline','grouporder',order')
set(gca,'YLim',[0,.6]);
ylabel('Allelic Fraction')




end


function test
maf='/Users/amaro/Downloads/Exomes.final_analysis_set.maf';
maf=load_table(maf);
 order={'024-CN-024-CN-N','006-CN-006-CN-N','003-CN-003-CN-N','008-CN-008-CN-N','001-CN-001-CN-N','009-CN-009-CN-N','005-CN-005-CN-N','022-CN-022-CN-N'...
,'025-CN-025-CN-N','026-CN-026-CN-N','027-CN-027-CN-N','028-CN-028-CN-N','002-CN-002-CN-N','029-CN-029-CN-N','007-CN-007-CN-N'}


Purity=alleleFreqPurity('/Users/amaro/Downloads/Exomes.final_analysis_purity.txt');

for i=1:size(order,2)
loc=find(ismember(maf.patient,order{i}),1);
sample_order{i,1}=maf.Tumor_Sample_Barcode{loc};
end

for i=1:length(sample_order)
loc=find(ismember(Purity.sample,sample_order{i}));
purity_values(i)=Purity.pur(loc);
purity_err(i)=Purity.purCI(loc,2);
end

afplots(maf,order,purity_values,purity_err)

end
