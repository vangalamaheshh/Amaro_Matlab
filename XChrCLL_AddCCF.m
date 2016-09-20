% Fill in X chr CCFs for CLL


%% Preprocessing 
gender=load_struct('/Users/amaro/Documents/CLL_Questions_For_Dan/XchrData_gender_assignment.txt');
P=load_table('/Users/amaro/Documents/StilgenbauerABSresults/Final_Stilgenbauer_2.4.ATW.ABSOLUTE.table.txt');
X0=load_table('/Users/amaro/Documents/StilgenbauerABSresults/FullStilgenbauerMafchrX.maf');
case_control_table=load_struct('/Users/amaro/Documents/StilgenbauerABSresults/case_control_tableStil.tsv');
AggregateAbsMaf=load_struct('/Users/amaro/Documents/StilgenbauerABSresults/AggregateStilgenbauerABS.maf');
AggregateAbsMaf.purity=str2double(AggregateAbsMaf.purity);
AggregateAbsMaf.ccf_hat=str2double(AggregateAbsMaf.ccf_hat);
AggregateAbsMaf.ccf_CI95_low=str2double(AggregateAbsMaf.ccf_CI95_low);
AggregateAbsMaf.ccf_CI95_high=str2double(AggregateAbsMaf.ccf_CI95_high);
X0=reorder_struct(X0,~isnan(X0.t_alt_count));
 [i m]=ismember(case_control_table.pair_id,P.array);
 P.case_sample(m(m>0),1)=case_control_table.case_sample(i);
% [i m]=ismember(gender.individual_id,X0.individual_id);
% X0.gender(m(m>0),1)=gender.gender(i);
f=fieldnames(AggregateAbsMaf);
ccf_bins={f{end-100:end}};
for i=1:slength(X0)
X0.individual_id{i,1}=X0.Tumor_Sample_Barcode{i}(1:13);
%X0.gender{i,1}=gender.gender{ismember(gender.individual_id,X0.individual_id{i})};
%X0.purity(i,1)=P.purity(ismember(P.case_sample,X0.Tumor_Sample_Barcode{i}));
end

X0=reorder_struct(X0,ismember(X0.Tumor_Sample_Barcode,P.case_sample));
%% Generate CCFs for Xchr mutations
CCF_range=0:0.01:2;

for i=1:slength(X0)
pur_X0=P.purity(ismember(P.case_sample,X0.Tumor_Sample_Barcode{i}));
gen_X0=gender.gender(ismember(gender.individual_id,X0.individual_id{i}));
if ~isempty(pur_X0)&&~isempty(gen_X0)
if strcmp(gen_X0,'Female')
CCFPdist=betapdf((CCF_range./2).*pur_X0,X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
else
CCFPdist=betapdf((CCF_range).*pur_X0,X0.t_alt_count(i)+1,X0.t_ref_count(i)+1);
end
CCFPdist=CCFPdist./sum(CCFPdist);
CCFPdist(101)=(sum(CCFPdist(102:end)))+CCFPdist(101);
[val,tmo]=max(CCFPdist);
total=sum(CCFPdist);

 X0.ccf_hat(i,1)=CCF_range(tmo);
 X0.purity(i,1)=pur_X0;
 X0.ccf_CI95_low(i,1)=CCF_range(find(cumsum(CCFPdist./sum(CCFPdist))>.05,1,'first'));
 X0.ccf_CI95_high(i,1)=CCF_range(find(cumsum(CCFPdist./sum(CCFPdist))<.95,1,'last'));
 for b=1:length(ccf_bins)
 X0.(ccf_bins{b}){i,1}=num2str(CCFPdist(b));
 end
else
    
   X0.CCF(i,1)=-1;
   X0.pci_high(i,1)=-1;
   X0.pci_low(i,1)=-1;
end
clear pur_X0
clear gen_X0
end

X0=rename_field(X0,'Tumor_Sample_Barcode','sample');


tmp=mergeStruct(AggregateAbsMaf,X0);
tmp=rmfield(tmp,'N');




