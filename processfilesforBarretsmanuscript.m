%process files for Barrett's paper

files=dir('~/Downloads/Sup_Table_2/');
cols_to_keep={'Hugo_Symbol';'Entrez_Gene_Id';'Center';'NCBI_Build';'Chromosome';'Start_position';'End_position';'Strand';'Variant_Classification';'Variant_Type';'Reference_Allele';'Tumor_Seq_Allele1';'Tumor_Seq_Allele2';'sample';'Protein_Change';'ref';'alt';'ccf_hat'};
for i=3:length(files)
    maf=load_struct(strcat(['~/Downloads/Sup_Table_2/',files(i).name]));
    for j=1:length(cols_to_keep)
        n_maf.(cols_to_keep{j})=maf.(cols_to_keep{j});
    end
    npieces=split(files(i).name,'-');
    numb=regexp(npieces{1},'UMBEER_','split');
    if isequal(npieces{2},'TP')
        name=strcat(numb{2},'-','Primary');
    else
        name=strcat(numb{2},'-','Barrett''s');
    end
    n_maf.sample=repmat({name},slength(n_maf),1);
    save_struct(n_maf,strcat('~/Downloads/Sup_Table_2_clean/',name,'.maf'));
end

files=dir('~/Downloads/Sup_Table_7');
ABS_Table=load_struct('~/Documents/ABSOLUTE_PP_SampleKey_clean.txt');
for i=5:length(files)
    maf=load_struct(strcat(['~/Downloads/Sup_Table_7/',files(i).name]));
    files(i).name
    for j=1:length(cols_to_keep)
        n_maf.(cols_to_keep{j})=maf.(cols_to_keep{j});
    end
    ntumorid=regexp(n_maf.sample{1},'Tumor-','split');
    k=~cellfun(@isempty,strfind(ABS_Table.sample,ntumorid{2}));
    n_maf.sample=repmat(ABS_Table.array(k),slength(n_maf),1);
    save_struct(n_maf,strcat('~/Downloads/Sup_Table_7_clean/',ABS_Table.array{k},'.maf'))
end