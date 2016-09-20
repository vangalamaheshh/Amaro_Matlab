sigABS=load_table('/Users/amaro/Documents/ProgressionMafUpdated3.16.maf');
for i=1:slength(sigABS)
    sigABS.indvidual{i,1}=sigABS.sample{i}(1:13);
end

sigABS.key=strcat(num2str(sigABS.Start_position),sigABS.indvidual);

genes=unique(sigABS.Hugo_Symbol);
for s=1:length(genes)
    all=[];
   g_maf=reorder_struct(sigABS,ismember(sigABS.Hugo_Symbol,genes{s})); 
   hold on   
    g_maf.ccf_CI95_high_hack=round(g_maf.ccf_CI95_high_hack*100)/100;

[n l]=count(sigABS.key(ismember(sigABS.Hugo_Symbol,genes{s})&ismember(sigABS.TP,1)));
    connect=l(n>1);
    all.key=l;
for i=1:slength(all)
    TPs=unique(g_maf.TP(ismember(g_maf.key,all.key{i})));
    all.ccf(i,1)=g_maf.ccf_mode_hack(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    all.ccf_high(i,1)=g_maf.ccf_CI95_high_hack(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    all.ccf_low(i,1)=g_maf.ccf_CI95_low_hack(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    if ismember(2,TPs)
    all.ccf2(i,1)=g_maf.ccf_mode_hack(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));
    all.ccf2_high(i,1)=g_maf.ccf_CI95_high_hack(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));
    all.ccf2_low(i,1)=g_maf.ccf_CI95_low_hack(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));

    else
        all.ccf2(i,1)=NaN;
        all.ccf2_high(i,1)=NaN;
        all.ccf2_low(i,1)=NaN;
    end
end

tp1x=.8:.4/slength(all):1.2;
tp2x=1.8:.4/slength(all):2.2;
all=sortstruct(all,'ccf',-1);
all.tp1x=[.8:.4/(slength(all)-1):1.2]';
all=sortstruct(all,'ccf2',-1);
all.tp2x=[1.8:.4/(slength(all)-1):2.2]';

all.ccf2_high(all.ccf2_high<all.ccf2)=all.ccf2(all.ccf2_high<all.ccf2);
all.ccf_high(all.ccf_high<all.ccf)=all.ccf(all.ccf_high<all.ccf);

 all.ccf_low=all.ccf-all.ccf_low;
 all.ccf2_low=all.ccf2-all.ccf2_low;
 all.ccf_high=all.ccf_high-all.ccf;
 all.ccf2_high=all.ccf2_high-all.ccf2;


for i=1:slength(all)
    

    TPs=unique(g_maf.TP(ismember(g_maf.key,all.key{i})));
    
   if ismember(2,TPs)&&ismember(1,TPs)
       if (all.ccf(i)+all.ccf_high(i))<(all.ccf2(i)-all.ccf2_low(i))
        errorbar([all.tp1x(i);all.tp2x(i)],[all.ccf(i);...
            all.ccf2(i)],[all.ccf_low(i);all.ccf2_low(i)],...
        [all.ccf_high(i);all.ccf2_high(i)],'r--')
       elseif (all.ccf(i)-all.ccf_low(i))>(all.ccf2(i)+all.ccf2_high(i))
                     errorbar([all.tp1x(i);all.tp2x(i)],[all.ccf(i);...
            all.ccf2(i)],[all.ccf_low(i);all.ccf2_low(i)],...
        [all.ccf_high(i);all.ccf2_high(i)],'b--')
       else 
 errorbar([all.tp1x(i);all.tp2x(i)],[all.ccf(i);...
            all.ccf2(i)],[all.ccf_low(i);all.ccf2_low(i)],...
        [all.ccf_high(i);all.ccf2_high(i)],'--','Color',[.5 .5 .5])
       end
               plot([all.tp1x(i)],all.ccf(i),'k.')
               plot([all.tp2x(i)],all.ccf2(i),'k.')

   elseif ismember(1,TPs)
      % plot([all.tp1x(i)],all.ccf(i),'k.')
   end


end
xlim([0.5 2.5])
ylim([0 1]);
title(strcat(genes{s}),'FontSize',28)
ylabel('Cancer Cell Fraction','FontSize',20)
set(gca,'YTick',[0;.5;1],'XTick',[1;2],'XTickLabel',{'TP1';'TP2'})
print(gcf,'-depsc',strcat('/Users/amaro/Documents/GeneProgressionPlotsCLLStilgenbauer/NoSolo/',genes{s},'.eps'));
close all

end