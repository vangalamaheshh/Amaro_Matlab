KIRC_maf=load_struct('~/Downloads/KIRC_aggregate_with_AF.maf');

code_class=get_coding_class_muts();
KIRC_maf.i_tumor_f=str2double(KIRC_maf.i_tumor_f);
KIRC_immune=load_struct('~/Documents/Rotations/Beck/KIRC_immune_cell_table.txt');

for i=1:slength(KIRC_data_matrix)
KIRC_data_matrix.nSigMuts(i,1)=sum(KIRC_maf.is_sig(ismember(KIRC_maf.tcga_id,KIRC_data_matrix.tcga_id{i})&ismember(KIRC_maf.Variant_Classification,code_class)));
KIRC_data_matrix.MATH(i,1)=100*mad(KIRC_maf.i_tumor_f(ismember(KIRC_maf.tcga_id,KIRC_data_matrix.tcga_id{i})&~isnan(KIRC_maf.i_tumor_f)))/median(KIRC_maf.i_tumor_f(ismember(KIRC_maf.tcga_id,KIRC_data_matrix.tcga_id{i})&~isnan(KIRC_maf.i_tumor_f)));
end
fs=fieldnames(KIRC_immune);
for i=1:slength(KIRC_data_matrix)
    for j=2:length(fs)
        KIRC_data_matrix.(fs{j})(i,1)=KIRC_immune.(fs{j})(ismember(KIRC_immune.PatientID,KIRC_data_matrix.tcga_id{i}));
    end
end

KIRC_seg_stats=load_struct('~/Downloads/SegFiles.tsv');
for i=1:slength(KIRC_seg_stats)
    KIRC_seg_stats.tcga_id{i,1}=KIRC_seg_stats.pair_id{i}(6:17);
end

for i=1:slength(KIRC_data_matrix)
    KIRC_data_matrix.varCNV{i,1}=KIRC_seg_stats.var_CNV{ismember(KIRC_seg_stats.tcga_id,KIRC_data_matrix.tcga_id{i})};
end
KIRC_CNV_matrix=load_struct('~/Documents/Rotations/Beck/CNV_matrix_KIRC.txt');
 for i=1:slength(KIRC_CNV_matrix)
KIRC_CNV_matrix.tcga_id{i,1}=KIRC_CNV_matrix.pair_id{i}(1:12);
 end
 fs=fieldnames(KIRC_CNV_matrix);
for i=1:slength(KIRC_data_matrix)
    for j=2:length(fs)-1
        if ~isempty(KIRC_CNV_matrix.(fs{j})(ismember(KIRC_CNV_matrix.tcga_id,KIRC_data_matrix.tcga_id{i})))
         KIRC_data_matrix.(fs{j})(i,1)=KIRC_CNV_matrix.(fs{j})(ismember(KIRC_CNV_matrix.tcga_id,KIRC_data_matrix.tcga_id{i}));
        else
        end
    end
end
KIRC_maf=reorder_struct(KIRC_maf ,ismember(KIRC_maf.Variant_Classification,code_class));
for i=1:slength(KIRC_data_matrix)
    for j=1:length(sig_genes)
        KIRC_data_matrix.(sig_genes{j})(i,1)=ismember(sig_genes{j},KIRC_maf.Hugo_Symbol(ismember(KIRC_maf.tcga_id,KIRC_data_matrix.tcga_id{i})));
    end
end

for i=1:slength(KIRC_data_matrix)
    if ~isempty(KIRC_clin.Status(ismember(KIRC_clin.tcga_id,KIRC_data_matrix.tcga_id{i})))
    KIRC_data_matrix.OS_Status(i,1)=KIRC_clin.Status(ismember(KIRC_clin.tcga_id,KIRC_data_matrix.tcga_id{i}));
    KIRC_data_matrix.days_surv(i,1)=KIRC_clin.Days(ismember(KIRC_clin.tcga_id,KIRC_data_matrix.tcga_id{i}));
    end
end

for i=1:slength(KIRC_data)
    if ~isempty(EST_KIRC.est_purity(ismember(EST_KIRC.barcode,KIRC_data.barcode{i})))
    KIRC_data.KIRC_Est_purity(i,1)=str2double(EST_KIRC.est_purity(ismember(EST_KIRC.barcode,KIRC_data.barcode{i})));
    end
end