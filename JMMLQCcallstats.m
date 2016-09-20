call_stats_files=load_struct('~/Projects/JMML2.0/NewSampleCallStatsQC/call_stats_files.txt');
outdir='~/Projects/JMML2.0/NewSampleCallStatsQC/';
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
    call.init_n_lod=str2double(call.init_n_lod);
    call.init_t_lod=str2double(call.init_t_lod);
    
    
    call=reorder_struct(call,~isnan(call.contig));
    all_calls=call;
    
    call=reorder_struct(call,((call.t_alt_count+call.t_ref_count)>30)&((call.n_alt_count+call.n_ref_count)>30));
     g = fittype('a*x')
    linfit=fit(call.tumor_f,call.normal_f,g)

    figure()
    hold on
    
    refline(linfit.a,0)
    plot(call.tumor_f,call.normal_f,'.','color',[190/255 190/255 190/255],'marker','o')
    ylim([0,1])
    x=call.tumor_f-call.normal_f;
     k=(call.tumor_f<.7)&(call.normal_f<(.40*call.tumor_f-.08));
    % k(ismember(call.judgement,'KEEP'))=1;
     plot(call.tumor_f(ismember(call.judgement,'KEEP')),call.normal_f(ismember(call.judgement,'KEEP')),'.','color',[50/255 190/255 50/255],'marker','o','LineWidth',2);
     sites_to_consider=(ismember(call.judgement,'KEEP')|...
         ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
          ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
          ismember(call.failure_reasons,'germline_risk'));
     
     if sum(k&sites_to_consider)>1
     plot(call.tumor_f(k&sites_to_consider),call.normal_f(k&sites_to_consider),'r.','MarkerSize',10)
     outfit=fit(call.tumor_f(k&sites_to_consider),call.normal_f(k&sites_to_consider),g);
     hline=refline(outfit.a,0);
     set(hline,'color','k')
     text(.1,.9,num2str(outfit.a));
     call=all_calls;
     
       for f=1:slength(call)
         call.logOdds(f,1)=normal_contamination_estimate_fitting(call.t_alt_count(f),call.t_ref_count(f)...
         ,call.n_alt_count(f),call.n_ref_count(f),outfit.a,linfit.a)  ;
       end
     
       normlod_filter=( call.logOdds>log(1000) & (ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
          ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
          ismember(call.failure_reasons,'germline_risk')));
                
                
      if sum(normlod_filter>0)
          
       call.judgement(normlod_filter)={'KEEP'};
       call_stats_files.TinN(i,1)=outfit.a;
       
       save_struct(call,strcat(outdir,call_stats_files.pair_id{i},'_call_stats_fixed.txt'));
      end
     end
     xlabel('Tumor Allele Fraction')
     ylabel('Normal Allele Fraction')
     title(call_stats_files.pair_id{i});
     saveas(gcf,strcat(outdir,call_stats_files.pair_id{i},'_NormalEst.fig')) 
     
     
   

     
end
