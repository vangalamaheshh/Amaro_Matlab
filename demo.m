  hold on   


[n l]=count(sigABS.key(ismember(sigABS.Hugo_Symbol,genes{s})&ismember(sigABS.TP,1)));
    connect=l(n>1);
    all.key=l;
for i=1:slength(all)
    TPs=unique(g_maf.TP(ismember(g_maf.key,all.key{i})));
    all.ccf(i,1)=g_maf.ccf_hat(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    if ismember(2,TPs)
    all.ccf2(i,1)=g_maf.ccf_hat(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));
    else
        all.ccf2(i,1)=NaN;
    end
end

tp1x=.8:.4/slength(all):1.2;
tp2x=1.8:.4/slength(all):2.2;
all=sortstruct(all,'ccf',-1);
all.tp1x=[.8:.4/(slength(all)-1):1.2]';
all=sortstruct(all,'ccf2',-1);
all.tp2x=[1.8:.4/(slength(all)-1):2.2]';

for i=1:slength(all)
    

    TPs=unique(g_maf.TP(ismember(g_maf.key,all.key{i})));
    
   if ismember(2,TPs)&&ismember(1,TPs)
       if all.ccf(i)<all.ccf2(i)
        plot([all.tp1x(i);all.tp2x(i)],[all.ccf(i);...
            all.ccf2(i)],'r--')
       elseif all.ccf(i)>all.ccf2(i)
             plot([all.tp1x(i);all.tp2x(i)],[all.ccf(i);...
            all.ccf2(i)],'b--')
       else 
           plot([all.tp1x(i);all.tp2x(i)],[all.ccf(i);...
            all.ccf2(i)],'--','Color',[.5 .5 .5])
       end
               plot([all.tp1x(i)],all.ccf(i),'k.')
               plot([all.tp2x(i)],all.ccf2(i),'k.')

   elseif ismember(1,TPs)
       plot([all.tp1x(i)],all.ccf(i),'k.')
   end


   
   
end


   xlim([0.5 2.5])

   title(strcat('Cancer Cell Fraction Shift for ',genes{s}),'FontSize',24)