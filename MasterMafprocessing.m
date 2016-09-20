MM=load_struct('~/Projects/CLL_Rush/MasterMaf/FINAL/AggregateMasterMaf_WithPoNFields.tsv');
SIF=load_struct('~/Projects/CLL_Rush/ABSOLUTE/CLL_filtered_group.txt');


MM.key=strcat(MM.sample,MM.Start_position,MM.Tumor_Seq_Allele2);
MM.validation_tumor_ref_count_targeted=str2double(MM.validation_tumor_ref_count_targeted);
MM.validation_tumor_alt_count_targeted=str2double(MM.validation_tumor_alt_count_targeted);
MM.comp_alt=ceil(MM.validation_tumor_alt_count_targeted/2);
MM.comp_ref=ceil(MM.validation_tumor_ref_count_targeted/2);
MM.discovery_tumor_alt_count_targeted=str2double(MM.discovery_tumor_alt_count_targeted);
MM.discovery_tumor_ref_count_targeted=str2double(MM.discovery_tumor_ref_count_targeted);
MM=reorder_struct(MM,~isnan(MM.validation_tumor_alt_count_targeted));
MM.discovery_tumor_alt_count_targeted(isnan(MM.discovery_tumor_alt_count_targeted))=0;
MM.discovery_tumor_ref_count_targeted(isnan(MM.discovery_tumor_ref_count_targeted))=0;

for i=8:slength(SIF)
    m=load_struct(SIF.maf_union_forcecalls{i});
    m.key=strcat(m.sample,m.Start_position,m.Tumor_Seq_Allele2);
    m.t_alt_count=str2double(m.t_alt_count);
    m.t_ref_count=str2double(m.t_ref_count);
    
    for j=1:slength(m)
        l=find(ismember(MM.key,m.key{j}));
        if ~isempty(l)
        if (MM.discovery_tumor_ref_count_targeted(l)+MM.discovery_tumor_alt_count_targeted(l))...
                >(MM.comp_ref(l)+MM.comp_alt(l))
            m.t_alt_count(j)=MM.discovery_tumor_alt_count_targeted(l);
            m.t_ref_count(j)=MM.discovery_tumor_ref_count_targeted(l);
        else
            m.t_alt_count(j)=MM.comp_alt(l);
            m.t_ref_count(j)=MM.comp_ref(l);
        end
        end
    end
   
    m=rmfield(m,'key');
    save_struct(m,SIF.maf_union_forcecalls{i});
end