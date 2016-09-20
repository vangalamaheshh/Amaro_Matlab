ICGC_RNA=load_table('/Users/amaro/Documents/ICGC_RNA_validation_results.maf');
min_val_count = 2;
ICGC_RNA.TotalTumorCov = ICGC_RNA.t_alt_count + ICGC_RNA.t_ref_count;
ICGC_RNA.TotalRNACov = ICGC_RNA.RNA_alt_count + ICGC_RNA.RNA_ref_count;

for i=1:slength(ICGC_RNA)
    ICGC_RNA.power_of_site(i,1)=1-hyge2cdf(min_val_count-1,ICGC_RNA.TotalRNACov(i),ICGC_RNA.t_alt_count(i),ICGC_RNA.TotalTumorCov(i));
end