SampleTable=load_struct('/xchip/cga-home/amaro/CLL/StilgenbauerRNAValidation/sample_id.rna_bam');
SampleTable=load_struct('/xchip/cga_home/amaro/CLL/StilgenbauerRNAValidation/sample_id.rna_bam');
SampleTable
SampleTable=reorder_struct(SampleTable,~ismember(SampleTable.bam_file_rna,''));
SampleTable
case_pair=load_struct('/xchip/cga_home/amaro/CLL/StilgenbauerRNAValidation/case_sample.pair_id.txt');

case_pair=reorder_struct(case_pair,ismember(case_pair.case_sample,SampleTable.sample_id));

pair_maf=load_struct('/xchip/cga_home/amaro/CLL/StilgenbauerRNAValidation/validated_mafs.txt');
pair_maf=reorder_struct(pair_maf,ismember(pair_maf.pair_id,case_pair.pair_id));
pair_maf
pair_maf=reorder_struct(pair_maf,~ismember(pair_maf.maf_file_capture_validated,''));
p=load_struct(pair_maf.maf_file_capture_validated{1});
P=p;
for i=2:slength(pair_maf)
p=load_struct(pair_maf.maf_file_capture_validated{i});
P=mergeStruct(p,P);
end

P=rmfield(P,'N');
sig_genes=load_struct('/xchip/cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/ICGC_StilgenbauerForComut1.16/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt');
for i=1:slength(pair_maf)
p=load_struct(pair_maf.indel_maf_file_capture_validated{i});
P=mergeStruct(p,P);
end
P=rmfield(P,'N');
P=reorder_struct(P,ismember(P.Hugo_Symbol,sig_genes.gene));
P.validation_judgement_rna=str2double(P.validation_judgement_rna);
P.validation_power_rna=str2double(P.validation_power_rna);