M=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1.maf');
[pon_vars pass_4 pass_classic]=get_pon_variant_type_with_allele_count_ll(M,'/cga/fh/pancan_data/pon/pon7/final_summed_tokens.hist.bin',3,-4);
M=reorder_struct(M,pass_4);

sigGeneMaf=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.4PoNCutOFF.sig_genes.maf');
sum(ismember(M.key,sigGeneMaf.key));
salvaged=reorder_struct(sigGeneMaf,~ismember(sigGeneMaf.key,M.key));
M=mergeStruct(M,salvaged);


case_control_table=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/case_control_table.tsv');
k=find(ismember(M.Matched_Norm_Sample_Barcode,case_control_table.case_sample));



for i=1:length(k)
loc=find(ismember(case_control_table.case_sample,M.Matched_Norm_Sample_Barcode{k(i)}));
M.Matched_Norm_Sample_Barcode{k(i)}=case_control_table.control_sample{loc};
end








% Just to confirm, were any SF3B1 mutations filtered due to PoN and not salvaged? 
M=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1.maf');
Ml=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.maf');

count(M.Protein_Change(~ismember(M.key,Ml.key)&ismember(M.Hugo_Symbol,'SF3B1')))
salvSF3=reorder_struct(M,~ismember(M.key,Ml.key)&ismember(M.Hugo_Symbol,'SF3B1')); 
Ml=mergeStruct(Ml,salvSF3);
 TP53Slav=reorder_struct(M,~ismember(M.key,Ml.key)&ismember(M.Hugo_Symbol,'TP53')); 
 Ml=mergeStruct(Ml,TP53Slav);
 Ml=rmfield(Ml,'N');
 save_struct(Ml,'~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.maf')

 
 IGV=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.4PoNCutOFF.sig_genes.txt');
 Ml=reorder_struct(Ml,~ismember(Ml.key,IGV.key(ismember(IGV.IGV0good1bad2notreviewed,'1'))));
 
 M=reorder_struct(M,~ismember(M.Tumor_Sample_Barcode,'ICGC_CLL-173-TD'));
 m178=load_struct('/local/cga-fh/cga/CLL_Wu_workspace/Pair/ICGC_CLL-178-TB-NBC--/jobs/capture/mut/oncotate/job.50063844/ICGC_CLL-178-TB-NBC--.snp.capture.maf.annotated')
 M=mergeStruct(M,m178);
 m178=load_struct('/local/cga-fh/cga/CLL_Wu_workspace/Pair/ICGC_CLL-178-TB-NBC--/jobs/m2/vcf2maf/job.50882290/ICGC_CLL-178-TD.maf');
m178.Matched_Norm_Sample_Barcode(:)={'ICGC_CLL-178-ND'};
m178.Matched_Norm_Sample_UUID(:)={'ICGC_CLL-178-ND'};
m178.Tumor_Sample_Barcode(:)={'ICGC_CLL-178-TD'};
M=mergeStruct(M,m178);

Mfull=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1.maf');
 