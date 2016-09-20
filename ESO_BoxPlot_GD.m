AggSeg=load_struct('/Users/amaro/Downloads/OLD_Eso_oncoprint/Agg_SegFile.txt');
Table=load_struct('/Users/amaro/Downloads/OLD_Eso_oncoprint/146_samples.ATW.ABSOLUTE.table.txt');

dS=reorder_struct(Table,~ismember(Table.Genomedoublings,'0'));
S=reorder_struct(Table,ismember(Table.Genomedoublings,'0'));

dSeg=reorder_struct(AggSeg,ismember(AggSeg.sample,dS.array));
ndSeg=reorder_struct(AggSeg,ismember(AggSeg.sample,S.array));
Table.Doubled(ismember(Table.Genomedoublings,'0'),1)={'Non-Doubled'};
Table.Doubled(~ismember(Table.Genomedoublings,'0'),1)={'Doubled'};
AggSeg.modal_total_cn=str2double(AggSeg.modal_total_cn);
Table.ploidy=str2double(Table.ploidy);
AggSeg.corrected_total_cn=str2double(AggSeg.corrected_total_cn);

for i=1:slength(Table)
    Table.HZ_segs(i,1)=sum(AggSeg.corrected_total_cn(ismember(AggSeg.sample,Table.array{i}))==0);
    
end


for i=1:slength(Table)
    Table.DeletedSegments(i,1)=sum(AggSeg.corrected_total_cn(ismember(AggSeg.sample,Table.array{i}))<=(mode(AggSeg.corrected_total_cn(ismember(AggSeg.sample,Table.array{i})))-1));
    Table.AmpSegs(i,1)=sum(AggSeg.corrected_total_cn(ismember(AggSeg.sample,Table.array{i}))./Table.ploidy(i)>2);
    Table.HAmpSegs(i,1)=sum(AggSeg.corrected_total_cn(ismember(AggSeg.sample,Table.array{i}))./Table.ploidy(i)>3);


end

distribution_ordered_plot(Table.DeletedSegments(~ismember(Table.Doubled,'Doubled')),Table.DeletedSegments(ismember(Table.Doubled,'Doubled')))
set(gca,'XTick',[1:2],'XTickLabel',{'Not Doubled';'Doubled'},'FontSize',18)
title('Number of Deleted Segments Doubled EAC versus Non-Doubled','FontSize',24)
ylabel('number of deleted segments','FontSize',18)
ranksum(Table.DeletedSegments(~ismember(Table.Doubled,'Doubled')),Table.DeletedSegments(ismember(Table.Doubled,'Doubled')))

distribution_ordered_plot(Table.HZ_segs(~ismember(Table.Doubled,'Doubled')),Table.HZ_segs(ismember(Table.Doubled,'Doubled')))
set(gca,'XTick',[1:2],'XTickLabel',{'Not Doubled';'Doubled'},'FontSize',18)
title('Number of HZ Deleted Segments Doubled EAC versus Non-Doubled','FontSize',24)
ylabel('number of homozygous deleted segments','FontSize',18)
ranksum(Table.HZ_segs(~ismember(Table.Doubled,'Doubled')),Table.HZ_segs(ismember(Table.Doubled,'Doubled')))

distribution_ordered_plot(Table.AmpSegs(~ismember(Table.Doubled,'Doubled')),Table.AmpSegs(ismember(Table.Doubled,'Doubled')))
set(gca,'XTick',[1:2],'XTickLabel',{'Not Doubled';'Doubled'},'FontSize',18)
title('Number of Amplified Segments Doubled EAC versus Non-Doubled','FontSize',24)
ylabel('number of amplfied segments','FontSize',18)
ranksum(Table.AmpSegs(~ismember(Table.Doubled,'Doubled')),Table.AmpSegs(ismember(Table.Doubled,'Doubled')))

distribution_ordered_plot(Table.HAmpSegs(~ismember(Table.Doubled,'Doubled')),Table.HAmpSegs(ismember(Table.Doubled,'Doubled')))
set(gca,'XTick',[1:2],'XTickLabel',{'Not Doubled';'Doubled'},'FontSize',18)
title('Number of High Amplified Segments Doubled EAC versus Non-Doubled','FontSize',24)
ylabel('number of highly amplfied segments','FontSize',18)
ranksum(Table.HAmpSegs(~ismember(Table.Doubled,'Doubled')),Table.HAmpSegs(ismember(Table.Doubled,'Doubled')))






SampleTable7=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/SevenMostRelatedPairs.txt');
ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Aggregate_Seg_file.txt');
ABS_seg_file.modal_total_cn=str2double(ABS_seg_file.modal_total_cn);
SampleTable7.ploidy=str2double(SampleTable7.ploidy);
ABS_seg_file.corrected_total_cn=str2double(ABS_seg_file.corrected_total_cn);

for i=1:slength(SampleTable7)
    SampleTable7.AmpSegs(i,1)=sum(((ABS_seg_file.corrected_total_cn(ismember(ABS_seg_file.sample,SampleTable7.pair_id{i})))./SampleTable7.ploidy(i))>3);
end
SampleTable7.tissue(~ismember(SampleTable7.tissue,'ESO'))={'Barretts'};

figure()
hold on
step=.2/(length(sort(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'ESO'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'ESO'))),25,e_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(mean(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'Barretts'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'Barretts'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(mean(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'Barretts'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

set(gca,'XTick',[2:3],'XTickLabel',{'Barretts';'EAC'},'FontSize',18)
ylabel('number of high level amplification segments','FontSize',18)
title('Number of high level amplifications in 7 most related pairs','FontSize',24)
ranksum(SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'Barretts')),SampleTable7.AmpSegs(ismember(SampleTable7.tissue,'ESO')))



ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Aggregate_Seg_file.txt');
SampleTable=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/BoxPlotWithDiagnosisSampleTable.txt');
absolute_table=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/BarrettsABSOLUTEtable.txt');

absolute_table.ploidy=str2double(absolute_table.ploidy);
ABS_seg_file.modal_total_cn=str2double(ABS_seg_file.modal_total_cn);

for i=1:slength(SampleTable)
    ploidy=absolute_table.ploidy(ismember(absolute_table.array,SampleTable.pair_id{i}));
    SampleTable.AmpSegs(i,1)=sum(((ABS_seg_file.corrected_total_cn(ismember(ABS_seg_file.sample,SampleTable.pair_id{i})))./ploidy)>2);
    SampleTable.DelSegs(i,1)=sum(ABS_seg_file.corrected_total_cn(ismember(ABS_seg_file.sample,SampleTable.pair_id{i}))<=(mode(ABS_seg_file.corrected_total_cn(ismember(ABS_seg_file.sample,SampleTable.pair_id{i})))-1));
end


    
    
e_clo=[34/255 139/255 34/255];
b_clo=[218/255 112/255 214/255];
ND=[255/255 182/255 193/255];

figure()
hold on
step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.tissue,'ESO'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.tissue,'ESO'))),25,e_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.tissue,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.tissue,'Dysplasia'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.tissue,'Dysplasia'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.tissue,'Dysplasia'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.tissue,'ND'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.tissue,'ND'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.tissue,'ND'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);
set(gca,'XTick',[1:3],'XTickLabel',{'ND';'Dysplasia';'EAC'},'FontSize',18)
ylabel('number of amplification segments','FontSize',18)
title('Number of amplifications in trios','FontSize',24)

figure()
hold on
step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.tissue,'ESO'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.tissue,'ESO'))),25,e_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.tissue,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.tissue,'Dysplasia'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.tissue,'Dysplasia'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.tissue,'Dysplasia'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.tissue,'ND'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.tissue,'ND'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.tissue,'ND'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);
set(gca,'XTick',[1:3],'XTickLabel',{'ND';'Dysplasia';'EAC'},'FontSize',18)
ylabel('number of deleted segments','FontSize',18)
title('Number of deletions in trios','FontSize',24)




figure()
hold on
step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO'))))-1);
scatter((4-.1):step:((4+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO'))),25,e_clo,'filled')
line([(4-.1):.05:(4+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'HGD'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'HGD'))),25,b_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'HGD'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'LGD'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'LGD'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'LGD'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'BE'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'BE'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'BE'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

set(gca,'XTick',[1:4],'XTickLabel',{'BE';'LGD';'HGD';'EAC'},'FontSize',18)
ylabel('number of amplification segments','FontSize',18)
title('Number of amplifications in trios','FontSize',24)

ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'HGD')))
ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'LGD')))
ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'BE')))

figure()
hold on
step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO'))))-1);
scatter((4-.1):step:((4+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO'))),25,e_clo,'filled')
line([(4-.1):.05:(4+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'HGD'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'HGD'))),25,b_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'HGD'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'LGD'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'LGD'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'LGD'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'BE'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'BE'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'BE'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

set(gca,'XTick',[1:4],'XTickLabel',{'BE';'LGD';'HGD';'EAC'},'FontSize',18)
ylabel('number of deleted segments','FontSize',18)
title('Number of deletions in trios','FontSize',24)

ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'HGD')))
ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'LGD')))
ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'BE')))




for i=1:slength(SampleTable)
    ploidy=absolute_table.ploidy(ismember(absolute_table.array,SampleTable.pair_id{i}));
    SampleTable.HighLevelAmpSegs(i,1)=sum((ABS_seg_file.corrected_total_cn(ismember(ABS_seg_file.sample,SampleTable.pair_id{i})))./ploidy>3);
    SampleTable.HZDelSegs(i,1)=sum(ABS_seg_file.corrected_total_cn(ismember(ABS_seg_file.sample,SampleTable.pair_id{i}))<=0);
end

figure()
hold on
step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ESO'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ESO'))),25,e_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'Dysplasia'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'Dysplasia'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'Dysplasia'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ND'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ND'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ND'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);
set(gca,'XTick',[1:3],'XTickLabel',{'ND';'Dysplasia';'EAC'},'FontSize',18)
ylabel('number of high level amplified segments','FontSize',18)
title('number of high level amplifications','FontSize',24)


figure()
hold on
step=.2/(length(sort(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ESO'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ESO'))),25,e_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'Dysplasia'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'Dysplasia'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'Dysplasia'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ND'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ND'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ND'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);
set(gca,'XTick',[1:3],'XTickLabel',{'ND';'Dysplasia';'EAC'},'FontSize',18)
ylabel('number of homozygous deletions','FontSize',18)
title('number of homozygous deletions in trios','FontSize',24)

ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'BE')))
ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'LGD')))
ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'HGD')))
ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'HGD')))
ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'LGD')))
ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'BE')))


ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.tissue,'ND')))
ranksum(SampleTable.DelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.DelSegs(ismember(SampleTable.tissue,'Dysplasia')))
ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.tissue,'ND')))
ranksum(SampleTable.AmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.AmpSegs(ismember(SampleTable.tissue,'Dysplasia')))

ranksum(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'ND')))
ranksum(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HighLevelAmpSegs(ismember(SampleTable.tissue,'Dysplasia')))

ranksum(SampleTable.HZDelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HZDelSegs(ismember(SampleTable.tissue,'ND')))
ranksum(SampleTable.HZDelSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HZDelSegs(ismember(SampleTable.tissue,'Dysplasia')))



boxplot(SampleTable.AmpSegs,SampleTable.tissue)
ylabel('number of amplified segments','FontSize',18)
title('Number of amplifications in trios','FontSize',24)
boxplot(SampleTable.DelSegs,SampleTable.tissue)
ylabel('number of deleted segments','FontSize',18)
title('Number of deletions in trios','FontSize',24)











figure()
hold on
step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO'))))-1);
scatter((4-.1):step:((4+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO'))),25,e_clo,'filled')
line([(4-.1):.05:(4+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'HGD'))))-1);
scatter((3-.1):step:((3+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'HGD'))),25,b_clo,'filled')
line([(3-.1):.05:(3+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'HGD'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'LGD'))))-1);
scatter((2-.1):step:((2+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'LGD'))),25,b_clo,'filled')
line([(2-.1):.05:(2+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'LGD'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

step=.2/(length(sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'BE'))))-1);
scatter((1-.1):step:((1+.1)),sort(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'BE'))),25,ND,'filled')
line([(1-.1):.05:(1+.1)],repmat(median(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'BE'))),5),'Color',[0 0 0],'LineStyle','--','LineWidth',2);

set(gca,'XTick',[1:4],'XTickLabel',{'BE';'LGD';'HGD';'EAC'},'FontSize',18)
ylabel('number of high level amplified segments','FontSize',18)
title('Number of high level amplifications in trios','FontSize',24)

ranksum(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'HGD')))
ranksum(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'LGD')))
ranksum(SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'ESO')),SampleTable.HighLevelAmpSegs(ismember(SampleTable.Diagnosis,'BE')))
