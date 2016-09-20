for i=1:slength(x)
maf=load_struct(x.col1{i});
trim_maf.position=maf.Start_position;
trim_maf.Chromosome=maf.Chromosome;
trim_maf.Hugo_Symbol=maf.Hugo_Symbol; 
trim_maf.homozygousix=maf.homozygousix;
trim_maf.sample=maf.sample;
trim_maf.clonalix=maf.clonalix;
trim_maf.q_hat=maf.q_hat;
trim_maf.modal_q_s=maf.modal_q_s;
trim_maf.HS_q_hat_1=maf.HS_q_hat_1;
trim_maf.HS_q_hat_2=maf.HS_q_hat_2;
trim_maf.ccf=maf.ccf_hat;
trim_maf.COSMIC_overlapping_mutations=maf.COSMIC_overlapping_mutations;
trim_maf.i_judgement=maf.i_judgement;
trim_maf.clonal_scna_mut_ix=maf.clonal_scna_mut_ix;
trim_maf.Variant_Classification=maf.Variant_Classification;
trim_maf.Variant_Type=maf.Variant_Type;
trim_maf.Tumor_Sample_UUID=maf.Tumor_Sample_UUID;

save_struct(trim_maf,strcat(x.col1{i},'_TRIM_MAF'));  
i
end
