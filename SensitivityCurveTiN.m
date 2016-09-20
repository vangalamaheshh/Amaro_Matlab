for ff=1:10

TruthData=load_struct('/Users/amaro/Downloads/STAD-TCGA-BR-8680-TP-NB-SM-3N2T7-SM-3N2TR.call_stats.txt');
TruthData.position=str2double(TruthData.position);


set(0,'DefaultFigureVisible','off')

[x, TinN(1,ff)]=TumorInNormalCallStats('/Users/amaro/Downloads/CallStatsTiN/8680_10_normal.callstats','8680_10percent_simulation',0,1.5,10^-4,100);

x=reorder_struct(x,ismember(x.position,TruthData.position));
TruthData=reorder_struct(TruthData,ismember(TruthData.position,x.position));
muts=ismember(TruthData.judgement,'KEEP');
k_1=ismember(x.judgement,'KEEP');
sens(1,ff)=sum(muts&k_1)/sum(muts);
spec(1,ff)=sum(~muts&~k_1)/sum(~muts);
N_mut(1,ff)=sum(muts);
TruthData_1=TruthData;
k=ismember(TruthData_1.judgement,'KEEP');

i=1;
while N_mut(i)>5
   i=i+1; 
l=find(k);
r_m=l(randi(length(l),5,1));
r_p=TruthData_1.position(r_m);



x=reorder_struct(x,~ismember(x.position,r_p));
TruthData_1=reorder_struct(TruthData_1,~ismember(TruthData_1.position,r_p));

k=ismember(TruthData_1.judgement,'KEEP');

N_mut(i,ff)=sum(k);
[x, TinN(i,ff)]=TumorInNormalCallStats(x,'8589_10percent_simulation',0,1.5,10^-4,100);
k_1=ismember(x.judgement,'KEEP');




slength(x)
slength(TruthData_1)
sens(i,ff)=sum(k&k_1)/sum(k);

spec(i,ff)=sum(~k&~k_1)/sum(~k);
end

end

counter=1;




for i=1:size(TinN,1)
    for j=1:size(TinN,2)
        

normlod_filter=( call.tumor_f>call.normal_f & call.logOdds>1.5 & call.diff_reads > .85 & (ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
        ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
        ismember(call.failure_reasons,'germline_risk')| ismember(call.failure_reasons,'alt_allele_in_normal,strand_artifact')));
    
    end
end