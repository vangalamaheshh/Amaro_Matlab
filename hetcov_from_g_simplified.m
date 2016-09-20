G_M=load_struct('~/Projects/Brain_Mets/Het_Cov/Germ_mafs.tsv');
G_M=reorder_struct(G_M,~ismember(G_M.germline_maf_singlesample_analysisready,''));
for i=1:slength(G_M)
    
    m=load_struct(G_M.germline_maf_singlesample_analysisready{i});
    m=reorder_struct(m,ismember(m.i_genotype,'0/1'));
    m=reorder_struct(m,ismember(m.Variant_Type,'SNP'));
    if mod(i,5)==0
        i
    end
    for j=1:slength(m)
        cov.Chromosome{j,1}=m.Chromosome{j,1};
        cov.Start_position{j,1}=m.Start_position{j,1};
        cov.Reference_Allele{j,1}=m.Reference_Allele{j,1};
        cov.Tumor_Seq_Allele1{j,1}=m.Tumor_Seq_Allele1{j,1};
        strs=split(m.i_allelic_depth{j},',');
        %cov.i_t_ref_count{j,1}=strs{1};
        %cov.i_t_alt_count{j,1}=strs{2};
        intervals.x1{j,1}=strcat(m.Chromosome{j},':',m.Start_position{j});
        
        
        
        
    end
    %save_struct(cov,strcat('/xchip/cga_home/amaro/Brain_Mets/Het_Cov/',G_M.pair_id{i},'.cov'));
    save_struct_noheader(intervals,strcat('/xchip/cga_home/amaro/Brain_Mets/Het_Cov/',G_M.pair_id{i},'.intervals'));
    G_M.intervals_hetsites{i,1}=strcat('/xchip/cga_home/amaro/Brain_Mets/Het_Cov/',G_M.pair_id{i},'.intervals');
    clear intervals
    
end