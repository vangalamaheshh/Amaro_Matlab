




deTiNResultsTable=load_table('~/Documents/deTiN_Figures/Figure2/PairsdeTiNInSilicoFigure2TablePanel1.tsv');
deTiNResultsTable=reorder_struct(deTiNResultsTable,~ismember(deTiNResultsTable.tumor_in_normal_allelic_CIs,''));


for i=1:slength(deTiNResultsTable)
 
    vals=regexp(deTiNResultsTable.tumor_in_normal_allelic_CIs{i},'  ','split');
    deTiNResultsTable.low_loh(i,1)=deTiNResultsTable.TiN_LoH(i)-str2double(vals{1});
    deTiNResultsTable.high_loh(i,1)=str2double(vals{2})-deTiNResultsTable.TiN_LoH(i);
    if ~isequal(deTiNResultsTable.tumor_in_normal_mut_ci{i},'')
    vals=regexp(deTiNResultsTable.tumor_in_normal_mut_ci{i},'  ','split');
    deTiNResultsTable.low_mut(i,1)=deTiNResultsTable.tumor_in_normal_estimate(i)-str2double(vals{1});
    deTiNResultsTable.high_mut(i,1)=str2double(vals{2})-deTiNResultsTable.tumor_in_normal_estimate(i);
     vals=regexp(deTiNResultsTable.combined_deTiN_CI_int{i},'  ','split');
    deTiNResultsTable.low_combined(i,1)=deTiNResultsTable.combined_deTiN_TiN_num(i)-str2double(vals{1});
    deTiNResultsTable.high_combined(i,1)=str2double(vals{2})-deTiNResultsTable.combined_deTiN_TiN_num(i);
    else
        deTiNResultsTable.low_mut(i,1)=0;
        deTiNResultsTable.high_mut(i,1)=0;
        deTiNResultsTable.low_combined(i,1)=deTiNResultsTable.TiN_LoH(i)-str2double(vals{1});
        deTiNResultsTable.high_combined(i,1)=str2double(vals{2})-deTiNResultsTable.TiN_LoH(i);
    end

    
        
end


figure()
hold on
errorbar(deTiNResultsTable.actual_TiN,deTiNResultsTable.combined_deTiN_TiN_num,deTiNResultsTable.low_combined,deTiNResultsTable.high_combined,'k.','MarkerSize',20,'Marker','square','LineStyle','none')

errorbar(deTiNResultsTable.actual_TiN,deTiNResultsTable.TiN_LoH,deTiNResultsTable.low_loh,deTiNResultsTable.high_loh,'b.','MarkerSize',15,'Marker','x','LineStyle','none')
errorbar(deTiNResultsTable.actual_TiN,deTiNResultsTable.tumor_in_normal_estimate,deTiNResultsTable.low_mut,deTiNResultsTable.high_mut,'r','MarkerSize',15,'Marker','x','LineStyle','none')

% plot(deTiNResultsTable.actual_TiN,deTiNResultsTable.TiN_LoH,'r.','MarkerSize',20)
% plot(deTiNResultsTable.actual_TiN,deTiNResultsTable.tumor_in_normal_estimate,'b.','MarkerSize',20)
xlim([0 1])
ylim([0 1])
legend('Combined','Allele Shift','Mutations')
xlabel('Actual Tumor In Normal','FontSize',20)
ylabel('Estimated Tumor In Normal','FontSize',20)





