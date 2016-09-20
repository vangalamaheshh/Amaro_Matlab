

SIF=load_struct('~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_4.22/AllelicCapsegOutput.tsv'); 


for i=1:slength(SIF)
seg_file=load_table(SIF.alleliccapseg_tsv{i}); 
seg_file=rmfield(seg_file,'header');
seg_file=rmfield(seg_file,'headline');
seg_file.sample=repmat(SIF.pair_id(i),slength(seg_file),1);
tumor_cov=load_table(SIF.Het_AD_Tumor{i});
tumor_cov=rmfield(tumor_cov,'header');
tumor_cov=rmfield(tumor_cov,'headline');
tumor_cov.Chromosome=chromosome2num_legacy(tumor_cov.Chromosome);
skew=0;
for j=1:slength(seg_file)
    hets=reorder_struct(tumor_cov,tumor_cov.Chromosome==seg_file.Chromosome(j)&tumor_cov.Start_position>=seg_file.Start_bp(j)&tumor_cov.Start_position<=seg_file.End_bp(j));
    if slength(hets)>1
    [seg_file.mu_minor(j,1), seg_file.mu_major(j,1), seg_file.sigma_major(j,1)]=allelic_sigma_mu_correction(hets,seg_file.tau(j));
    seg_file.sigma_minor(j,1)=seg_file.sigma_major(j,1);
    else
        seg_file.mu_major(j,1)=NaN;
        seg_file.mu_minor(j,1)=NaN;
        seg_file.sigma_major(j,1)=NaN;
        seg_file.sigma_minor(j,1)=NaN;
    end
        
end
save_struct(seg_file,strcat(['~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_4.22/',SIF.pair_id{i},'_mu_corrected.seg']))
SIF.seg_file_mu_fix{i,1}=strcat(['~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_4.22/',SIF.pair_id{i},'_mu_corrected.seg']);
i
end
save_struct(SIF,'~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_4.22/AllelicCapsegOutput.tsv')



SIF=load_struct('~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_5.17/ACS_sigma_SIF.txt');
ABS_table=load_table('~/Projects/CLL/ABSOLUTE/FullSetABSOLUTE.Table.CLL8.txt');
ABS_table=rmfield(ABS_table,'header');
ABS_table=rmfield(ABS_table,'headline');
for i=1:slength(SIF)
    ACS=load_table(SIF.ACS_sigma_fix_seg{i});
    ACS=rmfield(ACS,'header');
    ACS=rmfield(ACS,'headline');
    ACS=reorder_struct(ACS,~isnan(ACS.tau));
    k=find(ismember(ABS_table.array,ACS.sample{1}));
    X=SCNA_CCF_chisq(ACS,ABS_table.purity(k),ABS_table.ploidy(k),.005,1,2,.5,1,7);
    X=rmfield_if_exist(X,'CCF_prob');
    X=rmfield_if_exist(X,'CCF_prob_nextbest');
    save_struct(X,strcat(['~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_5.17/',SIF.pair_id{i},'_CCF.seg']))
    SIF.CCF_seg{i,1}=strcat(['~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_5.17/',SIF.pair_id{i},'_CCF.seg']);
    
    i
end

save_struct(SIF,'~/Projects/CLL/ABSOLUTE/SCNA_CCF_recall_5.17/AllelicCapsegOutput.tsv')
