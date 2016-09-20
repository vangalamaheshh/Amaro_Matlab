M=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1.maf');
ValidEvents=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ValidatedStilgenbauerEvents.txt');


external=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/external_id_capture');
[i m]=ismember(external.external_id_capture,ValidEvents.CollaboratorSampleID);
ValidEvents(m(m>0))=external.sample_id(i);



[pon_vars pass pass_classic]=get_pon_variant_type_with_allele_count_ll(M,'/cga/fh/pancan_data/pon/pon7/final_summed_tokens.hist.bin',3,-4);

SF3B1_mut_samples.TumorSample=ValidEvents.TumorName(ismember(ValidEvents.SF3B1MutationStatus1mutated,'1'));
count(pass(ismember(M.Tumor_Sample_Barcode,SF3B1_mut_samples.TumorSample)&ismember(M.Hugo_Symbol,'SF3B1')))
SF3B1_mut_samples.cDNA=ValidEvents.SF3B1cDNA(ismember(ValidEvents.SF3B1MutationStatus1mutated,'1'));
M.Protein_Change(~pass(ismember(M.Tumor_Sample_Barcode,SF3B1_mut_samples.TumorSample)&ismember(M.Hugo_Symbol,'SF3B1')))
M.Tumor_Sample_Barcode(~pass(ismember(M.Tumor_Sample_Barcode,SF3B1_mut_samples.TumorSample)&ismember(M.Hugo_Symbol,'SF3B1')))



TP53_mut_samples.TumorSample=ValidEvents.TumorName(ismember(ValidEvents.TP53Mutation1yes0no,'1'));
count(pass(ismember(M.Tumor_Sample_Barcode,TP53_mut_samples.TumorSample)&ismember(M.Hugo_Symbol,'TP53')))
TP53_mut_samples.cDNA=ValidEvents.TP53cDNA(ismember(ValidEvents.TP53Mutation1yes0no,'1'));
M.Protein_Change(~pass(ismember(M.Tumor_Sample_Barcode,TP53_mut_samples.TumorSample)&ismember(M.Hugo_Symbol,'TP53')))
M.Tumor_Sample_Barcode(~pass(ismember(M.Tumor_Sample_Barcode,TP53_mut_samples.TumorSample)&ismember(M.Hugo_Symbol,'TP53')))







pon_vars((ismember(M.Hugo_Symbol,'MYD88')&~pass&ismember(M.Protein_Change,'p.*160R')),:)

M=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1.maf');

[pon_vars pass_2_5 pass_classic]=get_pon_variant_type_with_allele_count_ll(M,'/cga/fh/pancan_data/pon/pon7/final_summed_tokens.hist.bin',3,-2.5);
[pon_vars pass_4 pass_classic]=get_pon_variant_type_with_allele_count_ll(M,'/cga/fh/pancan_data/pon/pon7/final_summed_tokens.hist.bin',3,-4);

count(M.Protein_Change(ismember(M.Hugo_Symbol,'SF3B1')&pass_2_5&~pass_4))
TP53Samples2=M.Tumor_Sample_Barcode(ismember(M.Hugo_Symbol,'TP53')&pass_2_5);
TP53Samples4=M.Tumor_Sample_Barcode(ismember(M.Hugo_Symbol,'TP53')&pass_4); 
TP53Samples2(~ismember(TP53Samples2,TP53Samples4))
% ans = 
%     'CLL-GCLL-0027-Tumor-SM-41JMQ'
%     'CLL-GCLL-0166-Tumor-SM-41JYX'
%     'CLL-GCLL-0202-Tumor-SM-4DP8D'
MIA=TP53Samples2(~ismember(TP53Samples2,TP53Samples4));


count(M.Protein_Change(ismember(M.Hugo_Symbol,'TP53')&pass_2_5&~pass_4&ismember(M.Tumor_Sample_Barcode,MIA)))
%       p.H193R: [1]
%       p.R158L: [1]
%       p.R175H: [1]












 sig=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.PoN4.0.MutSig/sig_genes4.0ll.txt');
 sig.q=str2double(sig.q);
 sig=reorder_struct(sig,sig.q<.1)
M=reorder_struct(M,ismember(M.Hugo_Symbol,sig.gene));
[pon_vars pass_4 pass_classic]=get_pon_variant_type_with_allele_count_ll(M,'/cga/fh/pancan_data/pon/pon7/final_summed_tokens.hist.bin',3,-4);
M.PoN_Status(pass_4,1)={'PASS'};
M.PoN_Status(~pass_4,1)={'SALVAGED'};

M=reorder_struct(M,ismember(M.Variant_Classification,'Nonsense_Mutation')|ismember(M.Variant_Classification,'Missense_Mutation')...
    |ismember(M.Variant_Classification,'Frame_Shift_Del')|ismember(M.Variant_Classification,'Frame_Shift_Ins')...
    |ismember(M.Variant_Classification,'In_Frame_Del')|ismember(M.Variant_Classification,'In_Frame_Ins')|...
    ismember(M.Variant_Classification,'Nonstop_Mutation')|ismember(M.Variant_Classification,'Splice_Site')...
    | ismember(M.Variant_Classification,'De_novo_Start_OutOfFrame')|ismember(M.Variant_Classification,'Start_Codon_SNP'));



[pon_vars pass_2 pass_classic]=get_pon_variant_type_with_allele_count_ll(M,'/cga/fh/pancan_data/pon/pon7/final_summed_tokens.hist.bin',3,-2); 
M=reorder_struct(M,pass_2);
