call_stats = load_struct('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan_NT_paper/recovery_sensitivity/aggregate_keep_call_stats_pancan.tsv');
tumors = unique(call_stats.tumor_name);
call_stats.xs=xhg19(chromosome2num_legacy(call_stats.contig),str2double(call_stats.position));
call_stats.tumor_f=str2double(call_stats.tumor_f);
pancan = load_table('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan_NT_paper/TiN_Estimates/estimates.tsv');
pancan = rmfield(pancan,{'header','headline'});  

for i = 1:slength(pancan)
pancan.normal_tissue{i,1}=pancan.pair_id{i}(end-1:end);
strs=regexp(pancan.pair_id{i},'_','split');
pancan.tumor_type{i,1}=strs{3};
pancan.tumor_sample{i,1}=strs{1};
pancan.normal_samples{i,1}=strs{2};
end

for i=1:length(tumors)
    blood_normal = pancan.normal_samples(ismember(pancan.tumor_sample,tumors{i})&ismember(pancan.normal_tissue,'NB'));
    figure_table.tumor{i,1} = tumors{i};
    figure_table.TiN_estimate(i,1) = pancan.combined_deTiN_TiN_num(ismember(pancan.tumor_sample,tumors{i})&ismember(pancan.normal_tissue,'NT'));
    sites = call_stats.xs(ismember(call_stats.normal_name,blood_normal)&call_stats.tumor_f>.1);
    figure_table.n_sites(i,1) = length(sites);
    deTiN_sites = call_stats.xs(ismember(call_stats.tumor_name,tumors{i})&~ismember(call_stats.normal_name,blood_normal));
    figure_table.deTiN_sensitivity(i,1) = sum(ismember(sites,deTiN_sites))/figure_table.n_sites(i);
    muTect_sites = call_stats.xs(ismember(call_stats.tumor_name,tumors{i})&~ismember(call_stats.normal_name,blood_normal)&ismember(call_stats.failure_reasons,''));
    figure_table.muTect_sensitivity(i,1) = sum(ismember(sites,muTect_sites))/figure_table.n_sites(i);
end

figure_table=reorder_struct(figure_table,~isnan(figure_table.TiN_estimate));
deTiN_vals = unique(figure_table.TiN_estimate);
deTiN_vals(deTiN_vals == .09) = [];
deTiN_vals(deTiN_vals == .91) = [];
for i = 1: length(deTiN_vals)
    
    aggregate_data(i,1) = mean(figure_table.muTect_sensitivity(figure_table.TiN_estimate==deTiN_vals(i)));
    aggregate_data(i,2) = mean(figure_table.deTiN_sensitivity(figure_table.TiN_estimate==deTiN_vals(i)));
    aggregate_data(i,3) = std(figure_table.muTect_sensitivity(figure_table.TiN_estimate==deTiN_vals(i)));
    aggregate_data(i,4) = std(figure_table.deTiN_sensitivity(figure_table.TiN_estimate==deTiN_vals(i)));
    
end

figure()
hold on
 plot(deTiN_vals,aggregate_data(:,2),'r--')
 plot(deTiN_vals,aggregate_data(:,1),'b--')
 errorbar(deTiN_vals,aggregate_data(:,1),aggregate_data(:,3),'b.','MarkerSize',20)
 errorbar(deTiN_vals,aggregate_data(:,2),aggregate_data(:,4),'r.','MarkerSize',20)
 xlim([0 1])
ylim([0 1])

%table averaging values to limit noise:
deTiN_summ_Table=[0,0.951883455000000,0.951883455000000,0.0537844530000000,0.0537844530000000;0.0100000000000000,0.923303335000000,0.928376511000000,0.0663604860000000,0.0675906020000000;0.0200000000000000,0.859970180000000,0.892155628000000,0.109014480000000,0.112466873000000;0.0300000000000000,0.719556170000000,0.858942648000000,0.112232797000000,0.113713042000000;0.0600000000000000,0.653808234000000,0.834750109000000,0.104619797000000,0.0653963840000000;0.136666667000000,0.355498962000000,0.730227787000000,0.0461973820000000,0.201858174000000;0.750000000000000,0.0869565220000000,0.173913043000000,0,0];
figure()
hold on
errorbar(deTiN_summ_Table(:,1),deTiN_summ_Table(:,2),deTiN_summ_Table(:,4),'b.','MarkerSize',20)
plot(deTiN_summ_Table(:,1),deTiN_summ_Table(:,2),'b--')
errorbar(deTiN_summ_Table(:,1),deTiN_summ_Table(:,3),deTiN_summ_Table(:,5),'r.','MarkerSize',20)
plot(deTiN_summ_Table(:,1),deTiN_summ_Table(:,3),'r--')
