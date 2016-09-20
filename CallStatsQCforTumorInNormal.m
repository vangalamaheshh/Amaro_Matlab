function CallStatsQCforTumorInNormal(call_stats_files,outdir)
call_stats_files=load_struct(call_stats_files);
for i=1:slength(call_stats_files)

    call=load_struct(call_stats_files.call_stats_capture{i});
    
    call.t_ref_count=str2double(call.t_ref_count);
    call.t_alt_count=str2double(call.t_alt_count);
    call.n_ref_count=str2double(call.n_ref_count);
    call.n_alt_count=str2double(call.n_alt_count);
    call.tumor_f=str2double(call.tumor_f);
    call.normal_f=str2double(call.normal_f);
    call.position=str2double(call.position);
    call.contig=chromosome2num(call.contig);
    call=reorder_struct(call,~isnan(call.contig));
    call=reorder_struct(call,((call.t_alt_count+call.t_ref_count)>50)&((call.n_alt_count+call.n_ref_count)>50));
     g = fittype('a*x')
    linfit=fit(call.tumor_f,call.normal_f,g)

    figure()
    hold on
    
    refline(linfit.a,0)
    plot(call.tumor_f,call.normal_f,'.','color',[190/255 190/255 190/255],'marker','o')
    ylim([0,1])
    x=call.tumor_f-call.normal_f;
     k=(call.tumor_f<.7)&(call.normal_f<(.40*call.tumor_f-.08));
     k(ismember(call.judgement,'KEEP'))=1;
     plot(call.tumor_f(ismember(call.judgement,'KEEP')),call.normal_f(ismember(call.judgement,'KEEP')),'.','color',[50/255 190/255 50/255],'marker','o','LineWidth',2);
     if sum(k)>1
         plot(call.tumor_f(k),call.normal_f(k),'r.','MarkerSize',10)
 outfit=fit(call.tumor_f(k),call.normal_f(k),g);
     hline=refline(outfit.a,0);
     set(hline,'color','k')
     text(.1,.9,num2str(outfit.a));
     end
     xlabel('Tumor Allele Fraction')
     ylabel('Normal Allele Fraction')
     title(call_stats_files.pair_id{i});
     saveas(gcf,strcat(outdir,call_stats_files.pair_id{i},'_NormalEst.fig')) 
end
end
function test

%call_stats_files=load_struct('~/Projects/JMML2.0/NewSampleCallStatsQC/call_stats_files.txt');


CallStatsQCforTumorInNormal('/Users/amaro/Downloads/CLLRushTinN.txt','/Users/amaro/Documents/CLL_RushQC')
end