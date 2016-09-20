function allelic_capseg_mu_sigma_correction_wrapper(seg_file,tumor_cov_file,pair_id)
seg_file=load_table(seg_file); 
seg_file=rmfields_if_exist(seg_file,'header');
seg_file=rmfields_if_exist(seg_file,'headline');
seg_file.sample=repmat({pair_id},slength(seg_file),1);


tumor_cov=load_table(tumor_cov_file);
tumor_cov=rmfields_if_exist(tumor_cov,'header');
tumor_cov=rmfields_if_exist(tumor_cov,'headline');
tumor_cov.Chromosome=chromosome2num_legacy(tumor_cov.Chromosome);


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

save_struct(seg_file,strcat([pair_id,'.corrected.sigma.seg']));

end