%% load MAFs

T1 = load_struct('/seq/hacohenlab1/dan/BTK/CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2_combinedmaf_dal.txt');
T2 = load_struct('/seq/hacohenlab1/dan/BTK/CLL-MDAC-0012-TS-NUNK-SM-3VIEE-SM-3VHZ2_combinedmaf.txt') ;

%% keep only coding

include_variants = { 'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Silent','Splice_Site'};

idx = zeros(slength(T1), 1);
for i = 1:length(include_variants)
    idx = idx | strcmp(T1.Variant_Classification, include_variants{i});
end
T1 = reorder_struct(T1, idx ==1);

clear vars idx 
idx = zeros(slength(T2), 1);
for i = 1:length(include_variants)
    idx = idx | strcmp(T2.Variant_Classification, include_variants{i});
end
T2 = reorder_struct(T2, idx ==1);

T1 = make_all_possible_numeric(T1);
T2 = make_all_possible_numeric(T2);

%% save if validated in either RNAseq or targeted seq or high AF 

idx1 = zeros(slength(T1),1);
idx1 = idx1 | T1.validation_tumor_alt_count_targeted>20;
idx1 = idx1 | T1.validation_tumor_alt_count_rna>5;
idx1 = idx1 | T1.i_tumor_f>0.15;

idx2 = zeros(slength(T2),1);
idx2 = idx2 | T2.validation_tumor_alt_count_targeted>20;
idx2 = idx2 | T2.validation_tumor_alt_count_rna>5;
idx2 = idx2 | T2.i_tumor_f>0.15;


%% salvage based on presence in the other sample

list_privilige_T1 = T1.Start_position(idx1);
for i = 1: length(list_privilige_T1)
    idx2 = idx2 | T2.Start_position==list_privilige_T1(i);
end

list_privilige_T2 = T2.Start_position(idx2);
for i = 1: length(list_privilige_T2)
    idx1 = idx1 | T1.Start_position==list_privilige_T2(i);
end



%% reject based on IGV
reject_list = {'ALG13', 'LZTS2','CHCHD10','NBPF10','WASH1','ZIC5','HLA-A', 'HLA-J', 'MST1P2', 'NBPF1','PRMT8','RPL13AP6','RPL32','RRN3', 'RBM33','RUNX1-IT1','WASH3P', 'MUC4','COL8A2', 'HES1'};


for i = 1:length(reject_list)
    idx = strcmp(T1.Hugo_Symbol, reject_list{i});
    idx1(idx,1) = 0;
    idx = strcmp(T2.Hugo_Symbol, reject_list{i});
    idx2(idx,1) = 0;
end
     

%% reject sites that were covered by validation effort and did not validated

idx = T1.validation_power_targeted>0.8 & T1.validation_judgement_targeted==0;

idx1(idx) =0;

idx = T2.validation_power_targeted>0.8 & T2.validation_judgement_targeted==0;

idx2(idx) =0;