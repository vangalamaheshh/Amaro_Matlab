simTable=load_table('~/Documents/bideTiNSimTable.txt');
errorbar(simTable.actual_mix,simTable.TiN_LoH,simTable.TiN_LoH-simTable.low_ci_loh,simTable.high_ci_loh-simTable.TiN_LoH,'r.')
hold on
errorbar(simTable.actual_mix,simTable.tumor_in_normal_estimate,simTable.tumor_in_normal_estimate-simTable.low_ci_mut,simTable.high_ci_mut-simTable.tumor_in_normal_estimate,'b.')

ylim([0 1]);
xlim([0 1])
h=refline(1,0);
set(h,'Color',[0 0 0],'LineStyle','--','LineWidth',1)

corrcoef(simTable.TiN_LoH,simTable.actual_mix)
corrcoef(simTable.tumor_in_normal_estimate(simTable.tumor_in_normal_estimate>=0),simTable.actual_mix(simTable.tumor_in_normal_estimate>=0))

figure()
color=jet(6);
inds=unique(simTable.individual_id);
hold on
ylim([0 1]);
xlim([0 1])
for i=1:length(inds)
errorbar(simTable.actual_mix(ismember(simTable.individual_id,inds{i}))...
    ,simTable.TiN_LoH(ismember(simTable.individual_id,inds{i})),...
    simTable.TiN_LoH(ismember(simTable.individual_id,inds{i}))-simTable.low_ci_loh(ismember(simTable.individual_id,inds{i}))...
    ,simTable.high_ci_loh(ismember(simTable.individual_id,inds{i}))-simTable.TiN_LoH(ismember(simTable.individual_id,inds{i})),'.','Color',color(i,:))

errorbar(simTable.actual_mix(ismember(simTable.individual_id,inds{i})),simTable.tumor_in_normal_estimate(ismember(simTable.individual_id,inds{i})),...
    simTable.tumor_in_normal_estimate(ismember(simTable.individual_id,inds{i}))-simTable.low_ci_mut(ismember(simTable.individual_id,inds{i})),...
    simTable.high_ci_mut(ismember(simTable.individual_id,inds{i}))-simTable.tumor_in_normal_estimate(ismember(simTable.individual_id,inds{i})),'.','Color',color(i,:))
end

h=refline(1,0);
set(h,'Color',[0 0 0],'LineStyle','--','LineWidth',1)

