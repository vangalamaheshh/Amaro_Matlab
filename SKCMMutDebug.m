five_SKCM=load_struct('~/Downloads/4_8_update_SKCM-TCGA-D3-A3ML-Tumor-SM-2SK4H_normal_mix_fiveSKCM-TCGA-D3-A3ML.snp.capture.maf.annotated');
fifty_SKCM=load_struct('~/Downloads/4_8_update_SKCM-TCGA-D3-A3ML-Tumor-SM-2SK4H_normal_mix_fiftySKCM-TCGA-D3-A3ML.snp.capture.maf.annotated');
pure=load_struct('~/Downloads/SKCM-TCGA-D3-A3ML-PureNormalSKCM-TCGA-D3-A3ML.snp.capture.maf.annotated');
pure=reorder_struct(pure,ismember(pure.Start_position,fifty_SKCM.Start_position));
pure=reorder_struct(pure,ismember(pure.Start_position,five_SKCM.Start_position));
five_SKCM=reorder_struct(five_SKCM,ismember(five_SKCM.Start_position,pure.Start_position));
fifty_SKCM=reorder_struct(fifty_SKCM,ismember(fifty_SKCM.Start_position,pure.Start_position));

fifty_SKCM.Start_position=str2double(fifty_SKCM.Start_position);
five_SKCM.Start_position=str2double(five_SKCM.Start_position);
pure.Start_position=str2double(pure.Start_position);

fifty_SKCM.i_normal_f=str2double(fifty_SKCM.i_normal_f);
five_SKCM.i_normal_f=str2double(five_SKCM.i_normal_f);
pure.i_normal_f=str2double(pure.i_normal_f);
fifty_SKCM.i_tumor_f=str2double(fifty_SKCM.i_tumor_f);
five_SKCM.i_tumor_f=str2double(five_SKCM.i_tumor_f);
pure.i_tumor_f=str2double(pure.i_tumor_f);

chr6q=64568368;
chr1p=120738866;

sum(ismember(fifty_SKCM.Chromosome,'1')&fifty_SKCM.Start_position<120738866)
sum(ismember(fifty_SKCM.Chromosome,'6')&fifty_SKCM.Start_position>64568368)

k=(ismember(fifty_SKCM.Chromosome,'1')&fifty_SKCM.Start_position<120738866);
j=(ismember(fifty_SKCM.Chromosome,'6')&fifty_SKCM.Start_position>64568368);
l=(ismember(fifty_SKCM.Chromosome,'7')|(ismember(fifty_SKCM.Chromosome,'9')&fifty_SKCM.Start_position>52844943));
boxplot([fifty_SKCM.i_normal_f(l)./fifty_SKCM.i_tumor_f(l);fifty_SKCM.i_normal_f(k)./fifty_SKCM.i_tumor_f(k);fifty_SKCM.i_normal_f(j)./fifty_SKCM.i_tumor_f(j)],[ones(sum(l),1);ones(sum(k),1).*2;ones(sum(j),1).*3])
