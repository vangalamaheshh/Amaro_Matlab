SIF=load_struct('~/Projects/CLL_Rush/ABSOLUTE/ABSOLUTE_mafs/DummySiF.txt');
m=load_struct(SIF.maf_union_forcecalls{1});

MM=m;

for i=2:slength(SIF)
    if isequal(SIF.Ind{i-1},SIF.Ind{i}) && i<slength(SIF)
        m=load_struct(SIF.maf_union_forcecalls{i});
        MM=mergeStruct(m,MM);
    else
        if i==slength(SIF)
            m=load_struct(SIF.maf_union_forcecalls{i});
            MM=mergeStruct(m,MM);
        end
         MM=rmfield(MM,'N');
         MM.t_alt_count=str2double(MM.t_alt_count);
        for j=1:slength(MM)
            k=ismember(MM.Start_position,MM.Start_position{j});
            if length(unique({MM.Tumor_Seq_Allele2{k}}))>1 && ~isequal(MM.Hugo_Symbol{j},'PLCG2') && length(unique({MM.i_judgement{k}}))>1
                MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'KEEP')&k));
                MM.Tumor_Seq_Allele1(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'KEEP')&k));
                MM.Variant_Classification(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Variant_Classification(ismember(MM.i_judgement,'KEEP')&k));
                MM.Protein_Change(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Protein_Change(ismember(MM.i_judgement,'KEEP')&k));
            elseif length(unique({MM.Tumor_Seq_Allele2{k}}))>1 && ~isequal(MM.Hugo_Symbol{j},'PLCG2')
                lls=find(k);
                [v p]=max(MM.t_alt_count(lls));
                
                MM.i_judgement{lls(p)}='KEEP';
                MM.i_failure_reasons{lls(p)}='';
                MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'KEEP')&k));
                MM.Tumor_Seq_Allele1(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'KEEP')&k));
                MM.Variant_Classification(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Variant_Classification(ismember(MM.i_judgement,'KEEP')&k));
                MM.Protein_Change(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Protein_Change(ismember(MM.i_judgement,'KEEP')&k));
            end
        end
       MM=reorder_struct(MM,~ismember(MM.Variant_Classification,'Silent'));
        samples=unique(MM.sample);
        for p=1:length(samples)
            m=reorder_struct(MM,ismember(MM.sample,samples{p}));
            
            save_struct(m,strcat('~/Projects/CLL_Rush/ABSOLUTE/ABSOLUTE_mafs/OutPutAlleleFix/',samples{p},'.tsv'));
        end
        if i<slength(SIF)
        clear MM
        m=load_struct(SIF.maf_union_forcecalls{i});
        MM=m;
        end
    end
end

%ANPEP