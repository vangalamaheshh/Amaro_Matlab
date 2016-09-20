 files = dir('/xchip/cga_home/amaro/CLL_Rush/ABSOLUTE/ABSOLUTE_mafs/');
files(1)=[];
files(1)=[];
MM=load_struct('/xchip/cga_home/amaro/CLL_Rush/ABSOLUTE/ABSOLUTE_mafs/MasterMaf_update.txt');

for i=1:length(files)
    m=load_struct(strcat('/xchip/cga_home/amaro/CLL_Rush/ABSOLUTE/ABSOLUTE_mafs/',files(i).name));
    m.key=strcat(m.sample,m.Chromosome,m.Start_position,m.Tumor_Seq_Allele2);
    for j=1:slength(m)
        
    k=ismember(MM.key,m.key{j});
    if isempty(find(k))
        m.i_judgement{j}='REJECT';
        m.i_failure_reasons{j}='fstar_tumor_lod';
    m.validation_tumor_ref_count_targeted{j,1}='0';
    m.validation_tumor_alt_count_targeted{j,1}='0';
    m.discovery_tumor_ref_count_targeted{j,1}='0';
    m.discovery_tumor_alt_count_targeted{j,1}='0';
    m.MM_judgement{j,1}='REJECT';
    m.MM_val_judgement{j,1}='0';
    else
    m.validation_tumor_ref_count_targeted{j,1}=MM.validation_tumor_ref_count_targeted{k};
    m.validation_tumor_alt_count_targeted{j,1}=MM.validation_tumor_alt_count_targeted{k};
    m.discovery_tumor_ref_count_targeted{j,1}=MM.discovery_tumor_ref_count_targeted{k};
    m.discovery_tumor_alt_count_targeted{j,1}=MM.discovery_tumor_alt_count_targeted{k};
    m.MM_judgement{j,1}=MM.i_judgement{k};
    m.MM_val_judgement{j,1}=MM.validation_judgement_targeted{k};
    
    if isequal(m.MM_val_judgement{j},'1')
        m.i_judgement{j}='KEEP';
        m.i_failure_reasons{j}='';
    elseif isequal(m.MM_val_judgement{j},'0')
        m.i_judgement{j}='REJECT';
        m.i_failure_reasons{j}='fstar_tumor_lod';
    elseif isequal(m.MM_val_judgement{j},'NaN')
%       cases  controls
%  pos   a        c       K
%  neg   b        d       -
% total  N        -       M

%  p = fexact( a,M,K,N, options)
%      M and N are constants described above. a and K are P-vectors. No checks for valid
%      input are made.
        M=str2double(m.validation_tumor_alt_count_targeted{j})+str2double( m.validation_tumor_ref_count_targeted{j})...
            +str2double(m.discovery_tumor_ref_count_targeted{j})+str2double(m.discovery_tumor_alt_count_targeted{j});
        K=str2double(m.validation_tumor_alt_count_targeted{j})+str2double(m.discovery_tumor_alt_count_targeted{j});
        N=str2double(m.validation_tumor_alt_count_targeted{j})+str2double( m.validation_tumor_ref_count_targeted{j});
        p_comp=fexact(str2double(m.validation_tumor_alt_count_targeted{j}),M,K,N)
        
        
        
        M=str2double(m.validation_tumor_alt_count_targeted{j})+str2double( m.validation_tumor_ref_count_targeted{j})...
            +str2double( m.validation_tumor_ref_count_targeted{j});
        K=str2double(m.validation_tumor_alt_count_targeted{j});
        N=str2double(m.validation_tumor_alt_count_targeted{j})+str2double( m.validation_tumor_ref_count_targeted{j});
        p_zero=fexact(str2double(m.validation_tumor_alt_count_targeted{j}),M,K,N)
        
        
        strcat(m.Hugo_Symbol{j},' ',m.Protein_Change{j})
        if  p_zero<.05
        m.i_judgement{j}='KEEP';
        m.i_failure_reasons{j}='';
        else
        m.i_judgement{j}='REJECT';
        m.i_failure_reasons{j}='fstar_tumor_lod';
        end
        
    end
    
    
    end

    end
    m=reorder_struct(m,~ismember(m.Hugo_Symbol,'GRN'));
    m=reorder_struct(m,~ismember(m.Hugo_Symbol,'SSTR1'));
    m=reorder_struct(m,~ismember(m.Hugo_Symbol,'GRIK2'));
    m=reorder_struct(m,~ismember(m.Hugo_Symbol,'SLC4A10'));
    m=reorder_struct(m,~ismember(m.Hugo_Symbol,'C19orf79'));
    m=reorder_struct(m,~ismember(m.Hugo_Symbol,'DCK'));
    
    

    save_struct(m,strcat('/xchip/cga_home/amaro/CLL_Rush/ABSOLUTE/ABSOLUTE_mafs/',files(i).name));

end


% for i=1:length(files)-1
%     if isequal({files(i).name(1:13)},{files(i+1).name(1:13)})
%     m=load_struct(strcat('/Users/amaro/Downloads/ABSOLUTE_mafs/11/',files(i+1).name));
%     MM=mergeStruct(MM,m);
%     else
%         pos=unique(MM.Start_position);
%         
%         for j=1:length(pos)
%             k=ismember(MM.Start_position,pos{j});
%             Alls=unique(MM.Tumor_Seq_Allele2(k));
%             if length(Alls)>1&&sum(k)<5
%             MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'KEEP')&k));
%             MM.Tumor_Seq_Allele1(ismember(MM.i_judgement,'REJECT')&k)=unique(MM.Tumor_Seq_Allele2(ismember(MM.i_judgement,'KEEP')&k));
%             
%             end
%         
%         end
%         
%         MM=rmfield(MM,'N');
%        if isfield(MM,'sample')
%        samples=unique(MM.sample);
%        for p=1:length(samples)
%            m=reorder_struct(MM,ismember(MM.sample,samples{p}));
%            save_struct(m,strcat('/Users/amaro/Downloads/ABSOLUTE_mafs/',samples{p},'.tsv'));
%        end
%        else 
%            pause=1;
%        end
%         
%         
%         
        
        
        
        
        
 