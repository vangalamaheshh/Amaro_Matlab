MM=load_struct('/Users/amaro/Downloads/MasterMafFixUp/MasterMaf.txt');

%files={'CLL-MDAC-0022-Tumor1-Normal.validated.maf';'CLL-MDAC-0022-Tumor2-Normal.validated.maf';'CLL-MDAC-0022-Tumor3-Normal.validated.maf';'CLL-MDAC-0022-Tumor4-Normal.validated.maf'};
files={'CLL-MDAC-0011-TP-NUNK-SM-3VI2P-SM-3VHZ1.validated.maf';
'CLL-MDAC-0011-TP-NUNK-SM-3VIEB-SM-3VHZ1.validated.maf';
'CLL-MDAC-0011-TP-NUNK-SM-4M92J-SM-3VHZ1.validated.maf';
'CLL-MDAC-0011-TP-NUNK-SM-4M92K-SM-3VHZ1.validated.maf';
'CLL-MDAC-0011-TS-NUNK-SM-3VIEC-SM-3VHZ1.validated.maf';
'CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2.validated.maf';
'CLL-MDAC-0012-TS-NUNK-SM-3VIEE-SM-3VHZ2.validated.maf';
};

MM.key=strcat(MM.sample,MM.Chromosome,MM.Start_position,MM.Tumor_Seq_Allele2);
%f=fieldnames(m);
%val_f={f{296:317}};
%val_f=val_f';
for i=1:length(files)
    m=load_struct(strcat('/Users/amaro/Downloads/MasterMafFixUp/',files{i,1}));
        m.key=strcat(m.sample,m.Chromosome,m.Start_position,m.Tumor_Seq_Allele2);
    for j=1:length(val_f)
        if isfield(m,val_f{j})
        for k=1:slength(m)
          
        l=find(ismember(MM.key,m.key{k}));
        if ~isempty(l)
            MM.(val_f{j}){l}=m.(val_f{j}){k};
        end
        end
        end
    end
    
    
end


files={'LN1-HS_ABS_MAF.txt';'L01-HS_ABS_MAF.txt';'BM1-HS_ABS_MAF.txt'};

for i=1:length(files)
    m=load_struct(strcat('/Users/amaro/Downloads/MasterMafFixUp/',files{i,1}));
        m.key=strcat(m.sample,m.Chromosome,m.Start_position);
        for k=1:slength(m)
          
        l=find(ismember(MM.key,m.key{k}));
        if ~isempty(l)
          MM.ccf_hat{l}=m.cancer_cell_frac{k};
        end
        end
end


MM=load_struct('/Users/amaro/Downloads/MasterMafFixUp/MasterMaf_update.txt');

MM.key=strcat(MM.sample,MM.Chromosome,MM.Start_position,MM.Tumor_Seq_Allele2);
files={'CLL-MDAC-0011-TP-NUNK-SM-3VIEB-SM-3VHZ1_ABS_MAF.txt';'CLL-MDAC-0011-TP-NUNK-SM-4M92J-SM-3VHZ1_ABS_MAF.txt';...
    'CLL-MDAC-0011-TP-NUNK-SM-4M92K-SM-3VHZ1_ABS_MAF.txt';'CLL-MDAC-0011-TS-NUNK-SM-3VIEC-SM-3VHZ1_ABS_MAF.txt';...
    'CLL-MDAC-0012-TP-NUNK-SM-3VIED-SM-3VHZ2_ABS_MAF.txt';'CLL-MDAC-0012-TS-NUNK-SM-3VIEE-SM-3VHZ2_ABS_MAF.txt';...
    'CLL-MDAC-0011-TP-NUNK-SM-3VI2P-SM-3VHZ1_ABS_MAF.txt'};

MM=reorder_struct(MM,~ismember(MM.sample,'NA'));
for i=1:length(files)
    m=load_struct(strcat('/Users/amaro/Downloads/MasterMafFixUp/',files{i,1}));
        m.key=strcat(m.sample,m.Chromosome,m.Start_position,m.Tumor_Seq_Allele2);
        fields=fieldnames(m);
        for k=1:slength(m)
          for f=1:length(fields)
              if isfield(MM,fields{f})
            l=find(ismember(MM.key,m.key{k}));
              
              if ~isempty(l)
                  MM.(fields{f}){l}=m.(fields{f}){k};
              end
              end
          end
        
          
        end
        
end


files=load_struct('~/Projects/CLL_Rush/ABSOLUTE/used_ABSOLUTE_mafs.txt');
for j=1:slength(files)
m=load_struct(files.maf_union_forcecalls{j});
m.t_alt_count=str2double(m.t_alt_count);                      
m.t_ref_count=str2double(m.t_ref_count);
m.t_f=(m.t_alt_count)./(m.t_alt_count+m.t_ref_count);
for i=1:slength(m)
if m.t_f(i)>.05 || m.t_alt_count(i)>9
m.i_failure_reasons{i}='';
m.i_judgement{i}='KEEP';
end
end
save_struct(m,strcat(files.maf_union_forcecalls{j},'.fix'));
end


