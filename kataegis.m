function kataegis(maf_file)

maf=load_struct(maf_file);

maf=reorder_struct(maf,ismember(maf.Variant_Type,'SNP'));

maf.x=xhg19(chrom2num(maf.Chromosome),str2double(maf.Start_position));

maf=sort_struct(maf,'x',1);

for i=1:slength(maf)
if isequal(maf.Reference_Allele{i},'T')||isequal(maf.Reference_Allele{i},'G')

    
    if isequal(maf.Reference_Allele{i},'T')
        maf.Reference_Allele{i}='A';

        if isequal(maf.Tumor_Seq_Allele2{i},'A')
        maf.Tumor_Seq_Allele2{i}='T';  
        elseif isequal(maf.Tumor_Seq_Allele2{i},'C')
        maf.Tumor_Seq_Allele2{i}='G'; 
        elseif isequal(maf.Tumor_Seq_Allele2{i},'G')
        maf.Tumor_Seq_Allele2{i}='C';
        end
    end

    if isequal(maf.Reference_Allele{i},'G')
        maf.Reference_Allele{i}='C';

        if isequal(maf.Tumor_Seq_Allele2{i},'A')
        maf.Tumor_Seq_Allele2{i}='T';  
        elseif isequal(maf.Tumor_Seq_Allele2{i},'C')
        maf.Tumor_Seq_Allele2{i}='G'; 
        elseif isequal(maf.Tumor_Seq_Allele2{i},'T')
        maf.Tumor_Seq_Allele2{i}='A';
        end
    end
end
end

maf.change=strcat(maf.Reference_Allele,'>',maf.Tumor_Seq_Allele2);
figure()


hold on
maf.distance=[diff(x);0];

    k=ismember(maf.change,'A>C');
    plot(x(k),maf.distance(k),'b.');
    
    k=ismember(maf.change,'A>T');
    plot(x(k),maf.distance(k),'.','color',[186/255 85/255 211/255]);
    
    k=ismember(maf.change,'A>G');
    plot(x(k),maf.distance(k),'g.');
    
    k=ismember(maf.change,'C>A');
    plot(x(k),maf.distance(k),'.','color',[135/255 206/255 250/255]);
     
    k=ismember(maf.change,'C>T');
    plot(x(k),maf.distance(k),'.','color',[1 215/255 0]);
    
    k=ismember(maf.change,'C>G');
    plot(x(k),maf.distance(k),'r.');
    set(gca,'yscale','log')



end

function test
maf_file='/Users/amaro/Downloads/STAD-TCGA-BR-7197-TP-NB-SM-2R9F5-SM-2R9FT.validated.maf';

kataegis(maf_file)
end
