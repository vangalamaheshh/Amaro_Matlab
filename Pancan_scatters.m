colheaders={'CRAN','MEN','MM','BCL','HN','SCLC'};
clear ylabel
scatter((.95+(1.05-.95).*rand(length(comparative_stats.Cranio_mutationrate),1)),comparative_stats.Cranio_mutationrate,35,[.35 .35 .35],'filled')
%set(gca,'YScale','log')
hold on
scatter((1.95+(2.05-1.95).*rand(length(comparative_stats.MEN_mutationrate),1)),comparative_stats.MEN_mutationrate',35,[.35 .35 .35],'filled')
scatter((2.95+(3.05-2.95).*rand(length(comparative_stats.MM_mutationrate),1)),comparative_stats.MM_mutationrate',35,[.35 .35 .35],'filled')
scatter((3.95+(4.05-3.95).*rand(length(comparative_stats.DLBCL_mutationrate),1)),comparative_stats.DLBCL_mutationrate,35,[.35 .35 .35],'filled')
scatter((4.95+(5.05-4.95).*rand(length(comparative_stats.HN_mutationrate),1)),comparative_stats.HN_mutationrate,35,[.35 .35 .35],'filled')
scatter((5.95+(6.05-5.95).*rand(length(comparative_stats.SCLC_mutationrate),1)),comparative_stats.SCLC_mutationrate,35,[.35 .35 .35],'filled')
set(gca,'XLim',[0.5,6.5])
set(gca,'XTick',[1,2,3,4,5,6])
set(gca,'XTickLabel',colheaders, 'FontSize',16);
ylabel('Nonsynonymous mutation rate/Mb','FontSize',18)
meds=[median(comparative_stats.Cranio_mutationrate),median(comparative_stats.MEN_mutationrate),median(comparative_stats.MM_mutationrate),median(comparative_stats.DLBCL_mutationrate),...
median(comparative_stats.HN_mutationrate),median(comparative_stats.SCLC_mutationrate)];
line([.8:.2:1.2],repmat(meds(1),3),'Color',[0,0,1])
line([1.8:.2:2.2],repmat(meds(2),3),'Color',[0,0,1])
line([2.8:.2:3.2],repmat(meds(3),3),'Color',[0,0,1])
line([3.8:.2:4.2],repmat(meds(4),3),'Color',[0,0,1])
line([4.8:.2:5.2],repmat(meds(5),3),'Color',[0,0,1])
line([5.8:.2:6.2],repmat(meds(6),3),'Color',[0,0,1])%  colheaders={'CRAN','MEN','MM','MEL'};
% 
% clear ylabel
% scatter((.95+(1.05-.95).*rand(length(comparative_stats.Mouse_dRanger),1)),comparative_stats.Mouse_dRanger',35,[.35 .35 .35],'filled')
% hold on
% %set(gca,'YScale','log')
% scatter((1.95+(2.05-1.95).*rand(length(comparative_stats.LUAD_dRanger),1)),comparative_stats.LUAD_dRanger,35,[.35 .35 .35],'filled')
% scatter((2.95+(3.05-2.95).*rand(length(comparative_stats.MM_dRanger),1)),comparative_stats.MM_dRanger,35,[.35 .35 .35],'filled')
% scatter((3.95+(4.05-3.95).*rand(length(comparative_stats.ME_dRanger),1)),comparative_stats.ME_dRanger,35,[.35 .35 .35],'filled')
% set(gca,'XLim',[0.5,4.5])
% set(gca,'XTick',[1,2,3,4])
% colheaders={'MSCLC','LUAD','MM','MEL'};
% set(gca,'XTickLabel',colheaders, 'FontSize',16);
% 
% ylabel('Number of rearrangements','FontSize',18)
% meds=[median(comparative_stats.Mouse_dRanger),median(comparative_stats.LUAD_dRanger)...
%     ,median(comparative_stats.MM_dRanger),median(comparative_stats.ME_dRanger)];
% line([.8:.2:1.2],repmat(meds(1),3),'Color',[0,0,1])
% line([1.8:.2:2.2],repmat(meds(2),3),'Color',[0,0,1])
% line([2.8:.2:3.2],repmat(meds(3),3),'Color',[0,0,1])
% line([3.8:.2:4.2],repmat(meds(4),3),'Color',[0,0,1])
% 


colheaders={'CRAN','MEN','MM','BCL'};

clear ylabel
scatter((.95+(1.05-.95).*rand(length(comparative_stats.Cranio_CN),1)),comparative_stats.Cranio_CN,35,[.35 .35 .35],'filled');
hold on
%set(gca,'YScale','log')
scatter((1.95+(2.05-1.95).*rand(length(comparative_stats.MEN_CN),1)),comparative_stats.MEN_CN,35,[.35 .35 .35],'filled');
scatter((2.95+(3.05-2.95).*rand(length(comparative_stats.MM_CN),1)),comparative_stats.MM_CN,35,[.35 .35 .35],'filled');
scatter((3.95+(4.05-3.95).*rand(length(comparative_stats.BCL_CN),1)),comparative_stats.BCL_CN,35,[.35 .35 .35],'filled');
set(gca,'XLim',[0.5,4.5])
set(gca,'XTick',[1,2,3,4])
set(gca,'XTickLabel',colheaders, 'FontSize',16);

ylabel('Fraction of genome affected by SCNAs','FontSize',18);


meds=[median(comparative_stats.Cranio_CN),median(comparative_stats.MEN_CN),...
    median(comparative_stats.MM_CN),median(comparative_stats.BCL_CN)];
line([.8:.2:1.2],repmat(meds(1),3),'Color',[0,0,1])
line([1.8:.2:2.2],repmat(meds(2),3),'Color',[0,0,1])
line([2.8:.2:3.2],repmat(meds(3),3),'Color',[0,0,1])
line([3.8:.2:4.2],repmat(meds(4),3),'Color',[0,0,1])
