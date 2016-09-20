

function [TiN TiNLod_Curve Somatic_Calls] =  Merged_SNV_LOH_TiN_estimation(call_stats_file,pair_id,firehose,iter_threshold,f_c,max_iter,chip,cyt_file,bed,cnv_blacklist)




[call TinN_mutect Curve]=TumorInNormalCallStats(call_stats_file,pair_id,firehose,iter_threshold,f_c,max_iter,chip);
call_pre=sum(ismember(call.failure_reasons,''));
[TiN,log_odds_curve,CI_l,CI_h]=TiNLOH(call_stats_file,'',firehose,cyt_file,bed,cnv_blacklist,'call_stats',pair_id);


MutC=Curve.Curve{end};
MutC=MutC-min(MutC);
mfix=MutC(1:2:201);
log_odds_curve=log_odds_curve-min(log_odds_curve);

combined=log_odds_curve+mfix';

n_log_c=exp(-combined)/sum(exp(-combined));
[CI_l_C,CI_h_C]=findCIp(n_log_c,.05,.001);
[v l]=min(combined);
TinN=l-1;
disp(sprintf('MuTect TiN %d CI: %d/%d',TinN_mutect,Curve.CI_l,Curve.CI_h));
disp(sprintf('Allelic TiN %d CI: %d/%d',TiN,CI_l,CI_h));

disp(sprintf('Combined TiN %d CI: %d/%d',TinN,CI_l_C,CI_h_C));

save(sprintf('TinN_%s.txt',pair_id),'TinN','-ascii', '-double')
save(sprintf('TinN_Mutect_%s.txt',pair_id),'TinN_mutect','-ascii', '-double')
save(sprintf('TinN_Allelic_%s.txt',pair_id),'TiN','-ascii', '-double')


 disp('calculating the log odds of being somatic of each possible mutation')
    for f=1:slength(call)
        if call.normal_f(f)>call.tumor_f(f)
        call.logOdds(f,1)=-Inf;
        call.prob(f,1)=0;
        else
        [call.logOdds(f,1) call.prob(f,1)]=normal_contamination_estimate_fitting(call.t_alt_count(f),call.t_ref_count(f)...
            ,call.n_alt_count(f),call.n_ref_count(f),TinN,1,.001)  ;
        end
    end

normlod_filter=( call.tumor_f>call.normal_f & call.logOdds>iter_threshold & call.diff_reads > .85 & (ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
        ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
        ismember(call.failure_reasons,'germline_risk')| ismember(call.failure_reasons,'alt_allele_in_normal')));
    
call.judgement(normlod_filter)={'KEEP'};


save_struct(call,strcat(pair_id,'.call_stats.TiN.txt'));
call_post=sum(ismember(call.judgement,'KEEP'));
save(sprintf('calls_pre_%s.txt',pair_id),'call_pre','-ascii', '-double');
save(sprintf('calls_post_%s.txt',pair_id),'call_post','-ascii', '-double');



figure()
hold on
plot(log_odds_curve,':','Color',[.5 .5 1])
plot(mfix,':','Color',[1 .5 .5])
plot(combined,'--','Color','k')

legend('Allelic','MuTecT','Combined')







end

function demo

%joint mutation and LOH TiN estimate demo
Merged_SNV_LOH_TiN_estimation('/Users/amaro/Downloads/8535_30N_70T.call_stats.pon_annotated.txt','8535-30-70','0'...
    ,.001,4,100,1.5,'/xchip/cga_home/amaro/TumorInNormal/FullGenomeCytoband.txt','NA',' ')
end