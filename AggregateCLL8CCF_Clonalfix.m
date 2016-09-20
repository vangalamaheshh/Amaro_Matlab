AggregateAbsMaf=load_struct('/Users/amaro/Documents/StilgenbauerABSresults/AggregateStilgenbauerABS.maf');
P=load_table('/Users/amaro/Documents/StilgenbauerABSresults/Final_Stilgenbauer_2.4.ATW.ABSOLUTE.table.txt');
case_control_table=load_struct('/Users/amaro/Documents/StilgenbauerABSresults/case_control_tableStil.tsv');

AggregateAbsMaf.q_hat=str2double(AggregateAbsMaf.q_hat);
AggregateAbsMaf.alt=str2double(AggregateAbsMaf.alt);
AggregateAbsMaf.ref=str2double(AggregateAbsMaf.ref);
[i m]=ismember(case_control_table.pair_id,P.array);
P.case_sample(m(m>0),1)=case_control_table.case_sample(i);
AggregateAbsMaf=reorder_struct(AggregateAbsMaf,ismember(AggregateAbsMaf.sample,P.case_sample));
CCF_range=0:0.01:2;
for i=1:slength(AggregateAbsMaf)
pur_X0=P.purity(ismember(P.case_sample,AggregateAbsMaf.sample{i}));
if AggregateAbsMaf.q_hat(i)>0
CCFPdist=betapdf((CCF_range./AggregateAbsMaf.q_hat(i)).*pur_X0,AggregateAbsMaf.alt(i)+1,AggregateAbsMaf.ref(i)+1);
CCFPdist=CCFPdist./sum(CCFPdist);
CCFPdist(101)=(sum(CCFPdist(102:end)))+CCFPdist(101);
[val,tmo]=max(CCFPdist);
 AggregateAbsMaf.ccf_new(i,1)=CCF_range(tmo);
 AggregateAbsMaf.purity_n(i,1)=pur_X0;
 AggregateAbsMaf.ccf_new_CI95_low(i,1)=CCF_range(find(cumsum(CCFPdist./sum(CCFPdist))>.05,1,'first'));
 AggregateAbsMaf.ccf_new_CI95_high(i,1)=max([CCF_range(find(cumsum(CCFPdist./sum(CCFPdist))<.95,1,'last')) 0]);
 
     
else
    CCFPdist=betapdf((CCF_range).*pur_X0,AggregateAbsMaf.alt(i)+1,AggregateAbsMaf.ref(i)+1);
    CCFPdist=CCFPdist./sum(CCFPdist);
    CCFPdist(101)=(sum(CCFPdist(102:end)))+CCFPdist(101);
    [val,tmo]=max(CCFPdist);
    AggregateAbsMaf.ccf_new(i,1)=CCF_range(tmo);
    AggregateAbsMaf.purity_n(i,1)=pur_X0;
    AggregateAbsMaf.ccf_new_CI95_low(i,1)=CCF_range(find(cumsum(CCFPdist./sum(CCFPdist))>.05,1,'first'));
    AggregateAbsMaf.ccf_new_CI95_high(i,1)=max([CCF_range(find(cumsum(CCFPdist./sum(CCFPdist))<.95,1,'last')) 0]);
end
if sum(CCFPdist(96:end))>.5
     AggregateAbsMaf.clonal_new(i,1)=1;
else
     AggregateAbsMaf.clonal_new(i,1)=0;

end
end