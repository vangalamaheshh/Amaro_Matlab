

Non_Cancer_RTK=sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND'))&(ismember(samples.EGFR,'1')|...
    ismember(samples.ERBB2,'1')|ismember(samples.ERBB3,'1')|ismember(samples.FGFR2,'1')|ismember(samples.VEGFA,'1')|...
    ismember(samples.KRAS,'1')|ismember(samples.ALK,'1')|ismember(samples.MET,'1')))/sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND')))


Cancer_RTK=sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO'))&(ismember(samples.EGFR,'1')|...
    ismember(samples.ERBB2,'1')|ismember(samples.ERBB3,'1')|ismember(samples.FGFR2,'1')|ismember(samples.VEGFA,'1')|...
    ismember(samples.KRAS,'1')|ismember(samples.ALK,'1')|ismember(samples.MET,'1')))/sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO')))


P=barh([Cancer_RTK,Non_Cancer_RTK]);
set(gca,'ycolor',[0 0 0],'ytick',[1,2],'yticklabel',{'Cancer';'Barretts'},'ylim',[0 3],'xlim',[0 1],'xtick',[0 0.5 1]);
xlabel('% altered'); box off
Pbaseline=get(P,'BaseLine');
set(Pbaseline,'Color',[.9 .9 .9],'LineWidth',1);
set(P(1),'facecolor',[255/255 127/255 80/255],'EdgeColor',[1,1,1])
title('RTK pathway')








Cancer_Cell=sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO'))&(...
    ismember(samples.CCND1,'1')|ismember(samples.CCNE1,'1')|ismember(samples.CDK6,'1')))/sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO')))

Non_Cancer_Cell=sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND'))&(...
    ismember(samples.CCND1,'1')|ismember(samples.CCNE1,'1')|ismember(samples.CDK6,'1')))/sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND')))

P=barh([Cancer_Cell,Non_Cancer_Cell]);
set(gca,'ycolor',[0 0 0],'ytick',[1,2],'yticklabel',{'Cancer';'Barretts'},'ylim',[0 3],'xlim',[0 1],'xtick',[0 0.5 1]);
xlabel('% altered'); box off
Pbaseline=get(P,'BaseLine');
set(Pbaseline,'Color',[.9 .9 .9],'LineWidth',1);
set(P(1),'facecolor',[255/255 127/255 80/255],'EdgeColor',[1,1,1])
title('Cell Cycle pathway')


Cancer_TF=sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO'))&(ismember(samples.MYB,'1')|...
    ismember(samples.GATA4,'1')|ismember(samples.GATA6,'1')|ismember(samples.MYC,'1')))/sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO')))

Non_Cancer_TF=sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND'))&(ismember(samples.MYB,'1')|...
    ismember(samples.GATA4,'1')|ismember(samples.GATA6,'1')|ismember(samples.MYC,'1')))/sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND')))


Non_Cancer_WNT=sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND'))&(ismember(samples.PIK3CA,'1')|...
    ismember(samples.CTNNB1,'1')))/sum((ismember(samples.tissue,'Dysplasia')|ismember(samples.tissue,'ND')))

Cancer_WNT=sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO'))&(ismember(samples.PIK3CA,'1')|...
    ismember(samples.CTNNB1,'1')))/sum((ismember(samples.tissue,'ESO')|ismember(samples.tissue,'ESO')))