sim_sif=load_table('~/Projects/TumorInNormal/SIF_for_deTiN_sims.txt');
sim_sif=rmfields_if_exist(sim_sif,'header');
sim_sif=rmfields_if_exist(sim_sif,'headline');

TiN_mix=[.995;.99;.98;.95;.93;.90;.80;.50;.25];
TiN_m_label={'half_p';'one';'two';'five';'seven';'ten';'twenty';'fifty';'seventy_five'};
for i=1:slength(sim_sif)
    for j=1:length(TiN_mix)
    unix(strcat(['bsub -q hour -M 10 -o stdout',num2str(i+j-1)...
        ,'.txt "java -jar ~/Projects/Tools/picard-tools-1.104/picard-tools-1.104/DownsampleSam.jar I='...
        ,sim_sif.bam_file_normal{i},' O=',sim_sif.sample_name{i},'_',num2str(TiN_mix(j)),'.bam'...
        ,' P=',num2str(TiN_mix(j)*sim_sif.ratio_coverage_tum_normal(i)),'"']));
    sim_sif.(strcat(['normal_mix_',TiN_m_label{j}])){i,1}=strcat([sim_sif.sample_name{i},'_',num2str(TiN_mix(j)),'.bam']);
    end
end
save_struct(sim_sif,'~/Projects/TumorInNormal/SIF_for_deTiN_sims.txt');


TiN_mix=[.005;.01;.02;.05;.07;.1;.2;.5;.75];
TiN_m_label={'half_p';'one';'two';'five';'seven';'ten';'twenty';'fifty';'seventy_five'};

for i=1:slength(sim_sif)
    for j=1:length(TiN_mix)
    unix(strcat(['bsub -q hour -M 10 -o stdout_tumor_mix',num2str(i+j-1)...
        ,'.txt "java -jar ~/Projects/Tools/picard-tools-1.104/picard-tools-1.104/DownsampleSam.jar I='...
        ,sim_sif.bam_file_tumor{i},' O=Tumor_mix_',sim_sif.sample_name{i},'_',num2str(TiN_mix(j)),'.bam'...
        ,' P=',num2str(TiN_mix(j)),'"']));
    sim_sif.(strcat(['tumor_mix_',TiN_m_label{j}])){i,1}=strcat(['Tumor_mix_',sim_sif.sample_name{i},'_',num2str(TiN_mix(j)),'.bam']);
    end
end

save_struct(sim_sif,'~/Projects/TumorInNormal/SIF_for_deTiN_sims.txt');


  n_fields={'normal_mix_half_p'
    'normal_mix_one'
    'normal_mix_two'
    'normal_mix_five'
    'normal_mix_seven'
    'normal_mix_ten'
    'normal_mix_twenty'
    'normal_mix_fifty'
    'normal_mix_seventy_five'
    };
    
    t_fields={'tumor_mix_half_p'
    'tumor_mix_one'
    'tumor_mix_two'
    'tumor_mix_five'
    'tumor_mix_seven'
    'tumor_mix_ten'
    'tumor_mix_twenty'
    'tumor_mix_fifty'
    'tumor_mix_seventy_five'};
    



for i=1:slength(sim_sif)
 for j=1:length(n_fields)
     
     unix(strcat(['bsub -M 10 -q week -o stdout_mixing_',num2str(i+j-1)...
         ,'.txt "java -jar ~/Projects/Tools/picard-tools-1.104/picard-tools-1.104/MergeSamFiles.jar I='...
         ,sim_sif.(n_fields{j}){i},' I=',sim_sif.(t_fields{j}){i},' O=MixedBam_',sim_sif.sample_name{i},'_',n_fields{j},'.bam"']))
    sim_sif.(strcat(['mixed_bam_','_',n_fields{j}])){i,1}=strcat(['MixedBam_',sim_sif.sample_name{i},'_',n_fields{j},'.bam']);
 end

end



mixed={'mixed_bam__normal_mix_half_p'
    'mixed_bam__normal_mix_one'
    'mixed_bam__normal_mix_two'
    'mixed_bam__normal_mix_five'
    'mixed_bam__normal_mix_seven'
    'mixed_bam__normal_mix_ten'
    'mixed_bam__normal_mix_twenty'
    'mixed_bam__normal_mix_fifty'
    'mixed_bam__normal_mix_seventy_five'};


for i=1:slength(sim_sif)
 for j=1:length(mixed)
     
     unix(strcat(['bsub -M 10 -q week -o stdout_addreplace',num2str(i+j-1),...
         ' "java -jar /xchip/cga_home/amaro/Tools/picard-tools-1.104/picard-tools-1.104/AddOrReplaceReadGroups.jar I=',...
         sim_sif.(mixed{j}){i},' O=ReadGroupFixed',sim_sif.(mixed{j}){i},' RGLB=Tumor RGPL=Illumina RGPU=barcode RGSM=',sim_sif.sample_name{i},'"']))
     
 end
end



