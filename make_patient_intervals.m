function make_patient_intervals(indel_maf,snp_maf,patient_ids)

%load firehose tables
patient_ids=load_struct(patient_ids);
indel_maf=load_struct(indel_maf);
snp_maf=load_struct(snp_maf);


%list of individuals
pats=unique(patient_ids.patient_id);


%get target annotation names
f=fieldnames(indel_maf);
f_ind=f{2};
f=fieldnames(snp_maf);
f_snp=f{2};

%remove samples missing SNP or indel mafs
disp('Removing the following pairs for missing SNVs:')
r_snv=snp_maf.pair_id(ismember(snp_maf.(f_snp),''));
r_snv

disp('Removing the following pairs for missing Indels:')
r_ind=indel_maf.pair_id(ismember(indel_maf.(f_ind),''));
r_ind

if ~isempty(r_ind)|~isempty(r_snv)
r_pairs=[r_snv;r_ind];
patient_ids=reorder_struct(patient_ids,~ismember(patient_ids.pair_id,r_pairs));
snp_maf=reorder_struct(snp_maf,~ismember(snp_maf.pair_id,r_pairs));
indel_maf=reorder_struct(indel_maf,~ismember(indel_maf.pair_id,r_pairs));
end

for i=1:length(pats)
    pat_index=find(ismember(patient_ids.patient_id,pats{i}));
    pats{i}
    %combine all SNPs and indels
    ind=load_struct(indel_maf.(f_ind){pat_index(1)});
    ind.key=strcat(ind.Chromosome,ind.Start_position,ind.Tumor_Seq_Allele2);
    snp=load_struct(snp_maf.(f_snp){pat_index(1)});
    snp.key=strcat(snp.Chromosome,snp.Start_position,snp.Tumor_Seq_Allele2);
    m=mergeStruct_ChecksIfEmpty(ind,snp);
    M=m;
    if length(pat_index)>1
    for j=2:length(pat_index)
        ind=load_struct(indel_maf.(f_ind){pat_index(j)});
        ind.key=strcat(ind.Chromosome,ind.Start_position,ind.Tumor_Seq_Allele2);
        snp=load_struct(snp_maf.(f_snp){pat_index(j)});
        snp.key=strcat(snp.Chromosome,snp.Start_position,snp.Tumor_Seq_Allele2);
        x=mergeStruct(snp,ind);
        if slength(x)<1
            xf=fieldnames(m);
            for f=1:length(xf)
            m(j).(xf{f}){1}=[];
            end
        else
        m(j)=mergeStruct_ChecksIfEmpty(ind,snp); 
        end
        M=mergeStruct_ChecksIfEmpty(M,m(j));
    end
    % M contains union of all calls in mafs from each pair
    % m contains a struct of structs for each pairs merged mafs
    M.key=strcat(M.Chromosome,M.Start_position,M.Tumor_Seq_Allele2);
    m=rmfield(m,'N');
    M=rmfield(M,'N');
   
    % generate aggregated maf file for each pair member of the patient
    % use keys to annotate forcecall status
    for j=1:length(pat_index)
    M_p=M;
    current_sample=unique(m(j).Tumor_Sample_Barcode);
    current_uuid=unique(m(j).Tumor_Sample_UUID);
    k=ismember(M_p.key,m(j).key);
    s=ismember(M_p.Tumor_Sample_Barcode,current_sample);
    
    % Scrub the maf so that it reflects just one pair
    M_p.i_failure_reasons=cell(slength(M_p),1);
    M_p.i_failure_reasons(~k)={'fstar_tumor_lod'};M_p.i_failure_reasons(k)={''};
    M_p.Tumor_Sample_Barcode(:)=current_sample;
    M_p.Tumor_Sample_UUID(:)=current_uuid;
    M_p.t_alt_count(~k)={'nan'};
    M_p.t_ref_count(~k)={'nan'};
    if isfield(M_p,'i_init_t_lod')
        M_p.i_init_t_lod(~k)={'nan'};
    end
    if isfield(M_p,'init_n_lod')
    M_p.init_n_lod(~k)={'nan'};
    end
    if isfield(M_p,'i_n_alt_count')
    M_p.i_n_alt_count(~k)={'nan'};
    end
    if isfield(M_p,'i_n_ref_count')
    M_p.i_n_ref_count(~k)={'nan'};
    end
    if isfield(M_p,'i_tumor_f')
    M_p.i_tumor_f(~k)={'nan'};
    end
    if isfield(M_p,'i_judgement')
    M_p.i_judgement(~k)={'REJECT'};
    end
    % Make interval list unique
    [num key]=count(M_p.key);
    dupes=key(find(num>1));
    dupelist=zeros(slength(M_p),1);
    for d=1:length(dupes)
        dupelist(find(ismember(M_p.key,dupes{d})&(~s)))=1;   
    end
    M_p=reorder_struct(M_p,~dupelist);
    
    save_struct(M_p,strcat(snp_maf.pair_id{pat_index(j)},'.union_patient_calls.maf'));
    end
    else
        % only one pair in this patient just rewrite the maf file with
        % failure reasons column added
        M.i_failure_reasons=cell(slength(M),1);M_p.i_failure_reasons(k)={''};
        save_struct(M,strcat(snp_maf.pair_id{pat_index(j)},'.union_patient_calls.maf'));
   
    end
    clear m M M_p
end
end
    
        
