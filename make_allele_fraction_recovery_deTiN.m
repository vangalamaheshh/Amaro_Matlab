


TableDeTiN=load_table('/Users/amaro/Documents/TiN_Table_inSiLiCO.txt');
mix=unique(TableDeTiN.mix);

SIF=load_struct('~/Documents/SIF_for_mutation_recovery_comparison.txt');
 individuals=unique(SIF.tumor);
 mixes=unique(SIF.mix);
 

 af=zeros(4,1);
 for i=1:length(individuals)
     pure_maf=load_struct(SIF.maf{ismember(SIF.tumor,individuals{i})&ismember(SIF.mix,mixes{1})});
     pure_maf.i_tumor_f=str2double(pure_maf.i_tumor_f);
     pure_maf.t_alt_count=str2double(pure_maf.t_alt_count);
     pure_maf.key=strcat(pure_maf.Start_position,pure_maf.Chromosome);
     pure_maf.t_ref_count=str2double(pure_maf.t_ref_count);
     
     pure_maf.total_reads=pure_maf.t_alt_count+pure_maf.t_ref_count;
     
     pure_maf=reorder_struct(pure_maf,(pure_maf.t_alt_count+pure_maf.t_ref_count)>20);
     SIF.nMutations(ismember(SIF.tumor,individuals{i})&ismember(SIF.mix,mixes{1}),1)=slength(pure_maf);
     af=quantile(pure_maf.i_tumor_f,3);
     pm_idx=find(ismember(SIF.tumor,individuals{i})&ismember(SIF.mix,mixes{1}));
     for j=1:length(mixes)
         mix_maf=load_struct(SIF.maf{ismember(SIF.tumor,individuals{i})&ismember(SIF.mix,mixes{j})});
         mm_indx=find(ismember(SIF.tumor,individuals{i})&ismember(SIF.mix,mixes{j}));                                       
         mix_maf.i_tumor_f=str2double(mix_maf.i_tumor_f); 
         mix_maf.t_alt_count=str2double(mix_maf.t_alt_count); 
         mix_maf.key=strcat(mix_maf.Start_position,mix_maf.Chromosome);
         data_matrix_total_default(i,j)=sum(ismember(mix_maf.key,pure_maf.key)&(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN')))/slength(pure_maf);
         data_matrix_total_deTiN(i,j)=sum(ismember(mix_maf.key,pure_maf.key))/slength(pure_maf);
         
         
         % do the same data matrix default and deTIN by different allele
         % fraction bins.
         nMuT_total(i,j)=slength(mix_maf);
         nMuT4(i,j)=sum(mix_maf.i_tumor_f>=af(3));
         nMuT3(i,j)=sum(mix_maf.i_tumor_f>=af(2)&mix_maf.i_tumor_f<=af(3));
         nMuT2(i,j)=sum(mix_maf.i_tumor_f>=af(1)&mix_maf.i_tumor_f<=af(2));
         nMuT1(i,j)=sum(mix_maf.i_tumor_f<=af(1));
         
         nMuT_totaldef(i,j)=sum((ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN')));
         nMuT4_def(i,j)=sum(mix_maf.i_tumor_f>=af(3)&(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN')));
         nMuT3_def(i,j)=sum(mix_maf.i_tumor_f>=af(2)&mix_maf.i_tumor_f<=af(3)&(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN')));
         nMuT2_def(i,j)=sum(mix_maf.i_tumor_f>=af(1)&mix_maf.i_tumor_f<=af(2)&(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN')));
         nMuT1_def(i,j)=sum(mix_maf.i_tumor_f<=af(1)&(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN')));
         
         nMuT4Pure(i,j)=sum(pure_maf.i_tumor_f>=af(3));
         nMuT3Pure(i,j)=sum(pure_maf.i_tumor_f>=af(2)&pure_maf.i_tumor_f<=af(3));
         nMuT2Pure(i,j)=sum(pure_maf.i_tumor_f>=af(1)&pure_maf.i_tumor_f<=af(2));
         nMuT1Pure(i,j)=sum(pure_maf.i_tumor_f<=af(1));
         
         
         data_matrix_4(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f>=af(3)),mix_maf.key))/sum(pure_maf.i_tumor_f>=af(3));
         data_matrix_3(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f>=af(2)&pure_maf.i_tumor_f<af(3)),mix_maf.key))/sum(pure_maf.i_tumor_f>=af(2)&pure_maf.i_tumor_f<af(3));
         data_matrix_2(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f>=af(1)&pure_maf.i_tumor_f<af(2)),mix_maf.key))/sum(pure_maf.i_tumor_f>=af(1)&pure_maf.i_tumor_f<af(2));
         data_matrix_1(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f<af(1)),mix_maf.key))/sum(pure_maf.i_tumor_f<af(1));
         
         data_matrix_4def(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f>=af(3)),mix_maf.key(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN'))))/sum(pure_maf.i_tumor_f>=af(3));
         data_matrix_3def(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f>=af(2)&pure_maf.i_tumor_f<af(3)),mix_maf.key(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN'))))/sum(pure_maf.i_tumor_f>=af(2)&pure_maf.i_tumor_f<af(3));
         data_matrix_2def(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f>=af(1)&pure_maf.i_tumor_f<af(2)),mix_maf.key(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN'))))/sum(pure_maf.i_tumor_f>=af(1)&pure_maf.i_tumor_f<af(2));
         data_matrix_1def(i,j)=sum(ismember(pure_maf.key(pure_maf.i_tumor_f<af(1)),mix_maf.key(ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN'))))/sum(pure_maf.i_tumor_f<af(1));
         
         
     end
 end
 hold on
 color=jet(8);
%  
% plot([1:10],mean(data_matrix_1(:,1:end)),'Color',color(1,:))
% plot([1:10],mean(data_matrix_2(:,1:end)),'Color',color(2,:))
% plot([1:10],mean(data_matrix_3(:,1:end)),'Color',color(3,:))
% plot([1:10],mean(data_matrix_4(:,1:end)),'Color',color(4,:))
% 
% plot([1:10],mean(data_matrix_1def(:,1:end)),'Color',color(5,:))
% plot([1:10],mean(data_matrix_2def(:,1:end)),'Color',color(6,:))
% plot([1:10],mean(data_matrix_3def(:,1:end)),'Color',color(7,:))
% plot([1:10],mean(data_matrix_4def(:,1:end)),'Color',color(8,:))

plot(str2double(mixes),median(data_matrix_total_default),'r--')
% errors

[phat pci_deTiN]=binofit(sum(nMuT4),max(sum(nMuT4)));
[phat pci_def]=binofit(sum(nMuT4_def),max(sum(nMuT4)));



[phat pci_deTiN_total]=binofit(sum(nMuT_total),max(sum(nMuT_total)));
[phat pci_def_total]=binofit(sum(nMuT_totaldef),max(sum(nMuT_totaldef)));

high_ci_4_deTiN=pci_deTiN(:,2)'-(sum(nMuT4)./max(sum(nMuT4)));
low_ci_4_deTin=(sum(nMuT4)./max(sum(nMuT4)))-pci_deTiN(:,1)';

high_ci_4_def=pci_def(:,2)'-(sum(nMuT4_def)./max(sum(nMuT4_def)));
low_ci_4_def=(sum(nMuT4_def)./max(sum(nMuT4_def)))-pci_def(:,1)';



high_ci_deTiN=pci_deTiN_total(:,2)'-(sum(nMuT_total)./max(sum(nMuT_total)));
low_ci_deTin=(sum(nMuT_total)./max(sum(nMuT_total)))-pci_deTiN_total(:,1)';

high_ci_def=pci_def_total(:,2)'-(sum(nMuT_totaldef)./max(sum(nMuT_totaldef)));
low_ci_def=(sum(nMuT_totaldef)./max(sum(nMuT_totaldef)))-pci_def_total(:,1)';

figure()
axes1 = axes('Parent',figure1,'XTick',[0.005 0.01 0.02 0.05 0.07 0.1 0.2],...
    'XScale','log',...
    'XMinorTick','on');
hold on
errorbar(str2double(mixes),(sum(nMuT4)./max(sum(nMuT4))),low_ci_4_deTin,high_ci_4_deTiN,'o-','Color',[215/255 25/255 28/255])
errorbar(str2double(mixes),(sum(nMuT4_def)./max(sum(nMuT4_def))),low_ci_4_def,high_ci_4_def,'o-','Color',[44/255 123/255 182/255])

errorbar(str2double(mixes),(sum(nMuT_total)./max(sum(nMuT_total))),low_ci_deTin,high_ci_deTiN,'o-','Color',[253/255 174/255 97/255])
errorbar(str2double(mixes),(sum(nMuT_totaldef)./max(sum(nMuT_totaldef))),low_ci_def,high_ci_def,'o-','Color',[171/255 217/255 233/255])


 xlim(axes1,[0 0.25]);
set(gca,'xscale','log')
hold on



% Change quartile 

xlabel('Tumor In Normal Mix','FontSize',25)
ylabel('Fraction Mutations Detected (N=??)','FontSize',25) %Range of mutations per sample
% driver figure? 

% add error bars 
legend({'DeTiN','No-DeTiN'}) % DeTiN versus No DeTiN


mix_maf=load_struct('~/Downloads/TCGA-BH-A0H7-TP-NT.snp.capture.maf.annotated');
mix_maf.x=xhg19(chromosome2num_legacy(mix_maf.Chromosome),str2double(mix_maf.Start_position));
mix_maf.i_normal_f=str2double(mix_maf.i_normal_f);
mix_maf.i_tumor_f=str2double(mix_maf.i_tumor_f);

mix_maf.i_n_alt_count=str2double(mix_maf.i_n_alt_count);
 mix_maf.i_n_ref_count=str2double(mix_maf.i_n_ref_count);

mix_maf.t_alt_count=str2double(mix_maf.t_alt_count);
mix_maf.t_ref_count=str2double(mix_maf.t_ref_count);
%mix_maf=reorder_struct(mix_maf,ismember(mix_maf.Chromosome,'17')|ismember(mix_maf.Chromosome,'16')|ismember(mix_maf.Chromosome,'18'));


k=ismember(mix_maf.i_failure_reasons,'')|ismember(mix_maf.i_failure_reasons,',PoN');

sites=find(~k);
for i=1:length(sites)
    
plot([mix_maf.x(sites(i)) mix_maf.x(sites(i))],[mix_maf.i_tumor_f(sites(i)) mix_maf.i_normal_f(sites(i))],'k--o','MarkerSize',10)
end
[phat pci]=binofit(mix_maf.t_alt_count,mix_maf.t_alt_count+mix_maf.t_ref_count,.32);
mix_maf.t_high_ci=pci(:,2)-mix_maf.i_tumor_f;
mix_maf.t_low_ci=mix_maf.i_tumor_f-pci(:,1);

[phat pci]=binofit(mix_maf.i_n_alt_count,mix_maf.i_n_alt_count+mix_maf.i_n_ref_count,.32);
mix_maf.n_high_ci=pci(:,2)-mix_maf.i_normal_f;
mix_maf.n_low_ci=mix_maf.i_normal_f-pci(:,1);
figure()
hold on
ha=errorbar(mix_maf.x,mix_maf.i_normal_f,mix_maf.n_low_ci,mix_maf.n_high_ci,'b.','MarkerSize',10,'LineWidth',1);
hb = get(ha,'children');  
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% Increase line length by 0.2 units
% Xdata(xleft) = Xdata(xleft) +40000000;
% Xdata(xright) = Xdata(xright) - 40000000;
% set(hb(2),'Xdata',Xdata)



ha=errorbar(mix_maf.x,mix_maf.i_tumor_f,mix_maf.t_low_ci,mix_maf.t_high_ci,'r.','MarkerSize',10,'LineWidth',1);
hb = get(ha,'children');  
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% Increase line length by 0.2 units
% Xdata(xleft) = Xdata(xleft) +40000000;
% Xdata(xright) = Xdata(xright) - 40000000;
% set(hb(2),'Xdata',Xdata)
% 
sites=find(~k);
% for i=1:length(sites)
% plot([mix_maf.x(sites(i)) mix_maf.x(sites(i))],[mix_maf.i_normal_f(sites(i)) mix_maf.i_tumor_f(sites(i))],'k--o','MarkerSize',10)
% end

plot(mix_maf.x(~k),mix_maf.i_tumor_f(~k),'r*','MarkerSize',10)
aL=num2chrom(16:19);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];

set(gca,'xtick',xM,'xticklabel',aL,'xlim',[min(xL)-10 max(xL)])
% tick lines to separate chromosomes 
line([xL xL]',[0+0*xL 1+0*xL]','color',0.25*[1 1 1],'linestyle','--','LineWidth',.5)
xlabel('Chromosome','FontSize',25)
ylabel('Allele fraction','FontSize',25)
legend({'Normal','Tumor','Recovered'})

figure()
m=load_struct(SIF.maf{1});
M=m;
for i=2:slength(SIF)
m=load_struct(SIF.maf{i});
M=mergeStruct(m,M);
end

k=ismember(M.i_failure_reasons,'')|ismember(M.i_failure_reasons,',PoN');
M.i_n_alt_count=str2double(M.i_n_alt_count);
[num alt_c]=count(M.i_n_alt_count(~k));

[num alt_c]=count(M.i_n_alt_count(k));
h=bar(str2double(alt_c),(num./sum(num)));
set(h,'FaceColor',[0 0 1])
pH = arrayfun(@(x) allchild(x),h);
set(pH,'FaceAlpha',0.5);

hold on
[num alt_c]=count(M.i_n_alt_count(~k));
h2=bar([1:10],sumns);
set(h2,'FaceColor',[1 0 0])
pH = arrayfun(@(x) allchild(x),h2);
set(pH,'FaceAlpha',0.5);

plot(str2double(alt_c),cumsum(num./sum(num)),'b--')
 