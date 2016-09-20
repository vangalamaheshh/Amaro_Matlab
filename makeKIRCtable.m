KIRC_purity=load_struct('~/Documents/SNP_array_map_pair_KIRC.purity.tsv');
KIRC_neoantigens=load_struct('~/Documents/Rotations/Beck/KIRC_neoantigens_total_mutations_etc.txt');
survival_data=load_struct('~/Documents/Rotations/Beck/KIRC_data/KIRC_clinical_core.txt');
OS=load_struct('~/Documents/Rotations/Beck/KIRC_data/KIRC_OS_core.txt');
for i=1:slength(KIRC_purity)
KIRC_purity.barcode{i,1}=KIRC_purity.pair_id{i}(6:17);
end
for i=1:slength(KIRC_neoantigens)
    x=ismember(KIRC_purity.barcode,KIRC_neoantigens.PatientID(i));
    KIRC_purity.neoantigen_observed(x,1)=KIRC_neoantigens.NeoAgs_ObservedExpected(i);
    KIRC_purity.neoantigen_expected(x,1)=KIRC_neoantigens.PredictedNeoAgs(i);
    KIRC_purity.total_mutations(x,1)=KIRC_neoantigens.TotalMutations(i);

end
% [i m]=ismember(sample_map.old,CNV_IGHV.samples);
% CNV_IGHV.samples(m(m>0))=sample_map.new(i);
KIRC_purity.grade=repmat({'NaN'},slength(KIRC_purity),1);
KIRC_purity.age=repmat({'NaN'},slength(KIRC_purity),1);
KIRC_purity.stage=repmat({'NaN'},slength(KIRC_purity),1);
KIRC_purity.gender=repmat({'NaN'},slength(KIRC_purity),1);
KIRC_purity.OS_y=repmat({'NaN'},slength(KIRC_purity),1);
KIRC_purity.OS_OS=repmat({'NaN'},slength(KIRC_purity),1);
for i=1:slength(survival_data)
    x=ismember(KIRC_purity.barcode,survival_data.feature(i));
    KIRC_purity.grade(x,1)=survival_data.grade(i);
    KIRC_purity.age(x,1)=survival_data.age(i);
    KIRC_purity.stage(x,1)=survival_data.stage(i);
    KIRC_purity.gender(x,1)=survival_data.gender(i);
end

for i=1:slength(OS)
    x=ismember(KIRC_purity.barcode,OS.feature(i));
    KIRC_purity.OS_y(x,1)=OS.OS_OS(i);
    KIRC_purity.OS_OS(x,1)=OS.OS_vital_status(i);
end

save_struct(KIRC_purity,'~/Documents/Rotations/Beck/MolecularDataKIRC.txt')