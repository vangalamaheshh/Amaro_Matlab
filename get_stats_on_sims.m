function [truepre post_calls actual_calls false_positives false_neg sen_pre sen_post]=get_stats_on_sims(TruthCallstats,MixedCallStatsPre,MixedCallStats)

T=load_struct(TruthCallstats);
M=load_struct(MixedCallStats);
Mre=load_struct(MixedCallStatsPre);


T=reorder_struct(T,ismember(T.position,M.position));
M=reorder_struct(M,ismember(M.position,T.position));
Mre=reorder_struct(Mre,ismember(Mre.position,T.position));

T=reorder_struct(T,ismember(T.contig,M.contig));
M=reorder_struct(M,ismember(M.contig,T.contig));
Mre=reorder_struct(Mre,ismember(Mre.contig,T.contig));

Mrep=Mre.position(ismember(Mre.judgement,'KEEP'));
Mp=M.position(ismember(M.judgement,'KEEP'));
Tp=T.position(ismember(T.judgement,'KEEP'));

pre_calls=length(Mrep);
post_calls=length(Mp);
actual_calls=length(Tp);
trueP=sum(ismember(Tp,Mp));
truepre=sum(ismember(Tp,Mrep));

false_positives=post_calls-trueP;
false_neg=actual_calls-trueP;
sen_pre=truepre/actual_calls;
sen_post=trueP/actual_calls;



end