call=load_struct('/Users/amaro/Downloads/CallStatsTiN/8589_10.call_stats.txt');
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
call.total_reads=str2double(call.total_reads);


TruthData=load_struct('/Users/amaro/Downloads/STAD-TCGA-BR-8589-TP-NB-SM-3N2T6-SM-3N2TQ.chr2.call_stats.txt');
TruthData.position=str2double(TruthData.position);
muts=ismember(TruthData.judgement,'KEEP');
true_muts=TruthData.position(muts);
contam_calls=ismember(call.judgement,'KEEP');
missed_positions=true_muts(~ismember(true_muts,call.position(contam_calls)));
f_ps=TruthData.position(ismember(TruthData.judgement,'REJECT'));
call=reorder_struct(call,ismember(call.position,TruthData.position));
TruthData=reorder_struct(TruthData,ismember(TruthData.position,call.position));


%storing total call stats

%limiting to 15x plus sites

 call.diff_reads=(call.t_alt_count+call.n_alt_count+call.t_ref_count+call.n_ref_count)./call.total_reads;
 %call.judgement(call.diff_reads <.85)={'REJECT'};
    counter=1;
for i=1:size(TinN,1)
    i
    for j=1:size(TinN,2)
        for f=1:slength(call)
        call.logOdds(f,1)=normal_contamination_estimate_fitting(call.t_alt_count(f),call.t_ref_count(f)...
            ,call.n_alt_count(f),call.n_ref_count(f),TinN(i,j),.99)  ;
        end

normlod_filter=( call.tumor_f>call.normal_f & call.logOdds>1.5  & (ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
        ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
        ismember(call.failure_reasons,'germline_risk')| ismember(call.failure_reasons,'alt_allele_in_normal,strand_artifact')));
    
    SEN_TiN(counter,1)=sum(ismember(missed_positions,call.position(normlod_filter)))/length(missed_positions);
    SPEC_TiN(counter,1)=sum(ismember(f_ps,call.position(~normlod_filter)))/sum(ismember(TruthData.judgement,'REJECT'));
    counter=counter+1;
    end
end