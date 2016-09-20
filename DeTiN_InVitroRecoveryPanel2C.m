
SIF=load_struct('~/Documents/SIF_for_In_Vitro_mutation_recovery.txt');
MForce=load_struct('~/Downloads/HCC-1143_100_N1.call_stats.txt');

individuals=unique(SIF.tumor);
mixes=unique(SIF.mix);

for j=1:length(individuals)
disp(sprintf('Analyzing Tumor %s',individuals{j}));
% MForce=load_struct(MSIF.fm{ismember(MSIF.tumor,individuals{j})});
MForce.n_alt_count=str2double(MForce.n_alt_count);
 MForce.n_ref_count=str2double(MForce.n_ref_count);
 MForce.normal_f=str2double(MForce.normal_f);
 
 pure_maf=load_struct(SIF.maf{ismember(SIF.mix,'0')&ismember(SIF.tumor,individuals{j})});
 %pure_maf=load_struct('/Volumes/xchip_cga_home/amaro/TumorInNormal/InVitroCRSP/mafs/HCC1143_calls.maf');
 %pure_maf=reorder_struct(pure_maf,ismember(pure_maf.Status,'TP'));
 pure_maf.key=strcat(pure_maf.Start_position,pure_maf.Chromosome);
 pure_maf.Start_position=str2double(pure_maf.Start_position);
% pure_maf=reorder_struct(pure_maf,ismember(pure_maf.Start_position,true_positives));
% zero_maf.key=strcat(zero_maf.Start_position,zero_maf.Chromosome);
% zero_maf=reorder_struct(zero_maf,ismember(zero_maf.key,pure_maf.key));
% pure_maf=zero_maf;
% pure_maf.t_alt_count=str2double(pure_maf.t_alt_count);
% pure_maf.t_ref_count=str2double(pure_maf.t_ref_count);
% pure_maf.total_reads=pure_maf.t_alt_count+pure_maf.t_ref_count;
% pure_maf=reorder_struct(pure_maf,pure_maf.total_reads>20);
%pure_maf.i_n_ref_count=str2double(pure_maf.i_n_ref_count);
pure_maf.i_tumor_f=str2double(pure_maf.i_tumor_f);
HighAF_Key=pure_maf.key(pure_maf.i_tumor_f>.3);
SIF.TotalMuts(ismember(SIF.tumor,individuals{j}),1)=slength(pure_maf);
SIF.TotalHighMuts(ismember(SIF.tumor,individuals{j}),1)=length(HighAF_Key);
TotalMuts=slength(pure_maf);
TotalHighMuts=length(HighAF_Key);
for mix=2:length(mixes)
    
    i=find(ismember(SIF.mix,mixes{mix})&ismember(SIF.tumor,individuals{j}));
    mix_maf=load_struct(SIF.maf{i});
    mix_maf.i_tumor_f=str2double(mix_maf.i_tumor_f);
    mix_maf.i_normal_f=str2double(mix_maf.i_normal_f);
    mix_maf.i_n_alt_count=str2double(mix_maf.i_n_alt_count);
    mix_maf.i_n_ref_count=str2double(mix_maf.i_n_ref_count);
    mix_maf.t_alt_count=str2double(mix_maf.t_alt_count);
    mix_maf.t_ref_count=str2double(mix_maf.t_ref_count);
    mix_maf=reorder_struct(mix_maf,~ismember(mix_maf.Chromosome,'NaN'));
    mix_maf.key=strcat(mix_maf.Start_position,mix_maf.Chromosome);
   % MForce = reorder_struct(MForce,~isnan(chromosome2num_legacy(MForce.contig)));
    MForce.key=strcat(MForce.position,MForce.contig);
    fpc_keys=mix_maf.key(~ismember(mix_maf.key,pure_maf.key));
    fpc_pos=zeros(length(fpc_keys),1);
   
    for z=1:length(fpc_keys)
        if mix_maf.i_n_alt_count(ismember(mix_maf.key,fpc_keys{z}))< MForce.n_alt_count(ismember(MForce.key,fpc_keys{z})) || MForce.normal_f(ismember(MForce.key,fpc_keys{z})) > mix_maf.i_normal_f(ismember(mix_maf.key,fpc_keys{z})) 
        fpc_pos(z)=0;
        else
            fpc_pos(z)=1;
        end
        
    end
    mix_maf=reorder_struct(mix_maf,~ismember(mix_maf.key,fpc_keys(fpc_pos==0)));
    SIF.muts_recovered(i,1)=sum(ismember(pure_maf.key,mix_maf.key));
    SIF.muts_recovered_no_deTiN(i,1)=sum(ismember(pure_maf.key,mix_maf.key(ismember(mix_maf.i_failure_reasons,{'',',PoN'}))));
    SIF.muts_recovered_high(i,1)=sum(ismember(pure_maf.key,mix_maf.key)&ismember(pure_maf.key,HighAF_Key));
    SIF.muts_recovered_high_no_deTiN(i,1)=sum(ismember(pure_maf.key,mix_maf.key(ismember(mix_maf.i_failure_reasons,{'',',PoN'})))&ismember(pure_maf.key,HighAF_Key));
    SIF.sensitivity(i,1)=sum(ismember(pure_maf.key,mix_maf.key))/TotalMuts;
    SIF.sensitivity_high_af(i,1)=sum(ismember(pure_maf.key,mix_maf.key)&ismember(pure_maf.key,HighAF_Key))/TotalHighMuts;

    SIF.false_positives(i,1)=sum(~ismember(mix_maf.key,pure_maf.key));
    SIF.false_positives_per_mutation(i,1)=sum(~ismember(mix_maf.key,pure_maf.key))/(SIF.muts_recovered(i,1)+SIF.false_positives(i,1));
    [l,SIF.lci(i,:)] =poissfit(SIF.false_positives(i,1),.32);
    
    [phat,SIF.pci(i,:)] = binofit(sum(ismember(pure_maf.key,mix_maf.key)),TotalMuts);
    [phat,SIF.pci_high(i,:)] = binofit(sum(ismember(pure_maf.key,mix_maf.key)&ismember(pure_maf.key,HighAF_Key)),TotalHighMuts);
    
end
end
SIF.mix=str2double(SIF.mix);
SIFp=reorder_struct(SIF,SIF.mix>0);

mixes=str2double(mixes)

Overall(1)=1;
Overall_high(1)=1;
NodeTiN(1)=1;
NodeTiN_high(1)=1;
PCI_overall(1,:)=[1,1];
PCI_overall_high(1,:)=[1,1];
PCI_NodeTiN(1,:)=[1,1];
PCI_NodeTiN_high(1,:)=[1,1];
for i=2:length(mixes)
    FDR_Overall(i-1)=median(SIF.false_positives_per_mutation(SIF.mix==mixes(i)));
    Overall(i)=sum(SIF.muts_recovered(SIF.mix==mixes(i)))/sum(SIF.TotalMuts(SIF.mix==mixes(i)));
    Overall_high(i)=sum(SIF.muts_recovered_high(SIF.mix==mixes(i)))/sum(SIF.TotalHighMuts(SIF.mix==mixes(i)));
    NodeTiN(i)=sum(SIF.muts_recovered_no_deTiN(SIF.mix==mixes(i)))/sum(SIF.TotalMuts(SIF.mix==mixes(i)));
    NodeTiN_high(i)=sum(SIF.muts_recovered_high_no_deTiN(SIF.mix==mixes(i)))/sum(SIF.TotalHighMuts(SIF.mix==mixes(i)));
    [phat,PCI_overall(i,:)]=binofit(sum(SIF.muts_recovered(SIF.mix==mixes(i))),sum(SIF.TotalMuts(SIF.mix==mixes(i))),.32);
    [phat,PCI_overall_high(i,:)]=binofit(sum(SIF.muts_recovered_high(SIF.mix==mixes(i))),sum(SIF.TotalHighMuts(SIF.mix==mixes(i))),.32);
    [phat,PCI_NodeTiN(i,:)]=binofit(sum(SIF.muts_recovered_no_deTiN(SIF.mix==mixes(i))),sum(SIF.TotalMuts(SIF.mix==mixes(i))),.32);
    [phat,PCI_NodeTiN_high(i,:)]=binofit(sum(SIF.muts_recovered_high_no_deTiN(SIF.mix==mixes(i))),sum(SIF.TotalHighMuts(SIF.mix==mixes(i))),.32);
    [phat,PCI_false_positives_per_mut(i,:)]=binofit(sum(SIF.false_positives(SIF.mix==mixes(i))),sum(SIF.TotalMuts(SIF.mix==mixes(i))),.32);

    
end
    figure()
       SIFp=sort_struct(SIFp,'mix')
    errorbar(mixes(2:end),FDR_Overall,FDR_Overall(1:end)-PCI_false_positives_per_mut(2:end,1)',PCI_false_positives_per_mut(2:end,2)'-FDR_Overall(1:end),'k.','LineStyle','-','MarkerSize',30,'Color',[.5 .5 .5])
    %errorbar(SIFp.mix,SIFp.false_positives_per_mutation,SIFp.false_positives_per_mutation-PCI_false_positives_per_mut(2:end,1),PCI_false_positives_per_mut(2:end,2)-(SIFp.false_positives_per_mutation),'k.','LineStyle','--','MarkerSize',30)
hold on
    errorbar(mixes(1:end),Overall(1:end),Overall(1:end)-PCI_overall(1:end,1)',PCI_overall(1:end,2)'-Overall(1:end),'.','LineStyle','-','MarkerSize',30,'Color',[202/255,0,32/255])
    
    errorbar(mixes(1:end),NodeTiN(1:end),NodeTiN(1:end)-PCI_NodeTiN(1:end,1)',PCI_NodeTiN(1:end,2)'-NodeTiN(1:end),'.','LineStyle','-','MarkerSize',30,'Color',[5/255,113/255,176/255])
    xlim([0 1])
    ylim([0 1])
ylabel('False Positive Rate','FontSize',30)
xlabel('Tumor In Normal Fraction','FontSize',30)

figure()
errorbar(mixes(2:end),Overall(2:end),Overall(2:end)-PCI_overall(2:end,1)',PCI_overall(2:end,2)'-Overall(2:end),'.','LineStyle','None','MarkerSize',30,'Color',[202/255,0,32/255])
hold on
errorbar(mixes(2:end),Overall_high(2:end),Overall_high(2:end)-PCI_overall_high(2:end,1)',PCI_overall_high(2:end,2)'-Overall_high(2:end),'.','LineStyle','-','MarkerSize',30,'Color',[244/255,165/255,130/255])
errorbar(mixes(2:end),NodeTiN(2:end),NodeTiN(2:end)-PCI_NodeTiN(2:end,1)',PCI_NodeTiN(2:end,2)'-NodeTiN(2:end),'.','LineStyle','None','MarkerSize',30,'Color',[5/255,113/255,176/255])
errorbar(mixes(2:end),NodeTiN_high(2:end),NodeTiN_high(2:end)-PCI_NodeTiN_high(2:end,1)',PCI_NodeTiN_high(2:end,2)'-NodeTiN_high(2:end),'.','LineStyle','-','MarkerSize',30,'Color',[146/255,197/255,222/255])

errorbar(SIFp.mix,SIFp.sensitivity_high_af,SIFp.sensitivity_high_af-SIFp.pci_high(:,1),SIFp.pci_high(:,2)-SIFp.sensitivity_high_af,'r.','LineStyle','-','MarkerSize',30)
legend('all mutations','mut > 30% af')
xlabel('Tumor In Normal Fraction','FontSize',30)
ylabel('Fraction of mutations recovered','FontSize',30)
