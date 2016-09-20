% gene_exons=load_struct('/Users/amaro/Downloads/gencode_basic_exons.merged_new.txt');
% gene_exons.start=str2double(gene_exons.start);
% gene_exons.end=str2double(gene_exons.end);
% 
% 
% for i=1:slength(gene_exons)
% str=split(gene_exons.gene_name{i},';');
% gene_exons.name2{i,1}=str{1};
% end


start_g=1;
g_count=1;
counter=1;


while counter <=slength(gene_exons)
    c_gene=gene_exons.name2{counter};
    c_chr=gene_exons.cur{counter};
    
    if counter < slength(gene_exons)
        if isequal(c_gene,gene_exons.name2{counter+1})
            counter=counter+1;
        else
            gene_merge.gene{g_count,1}=c_gene;
            gene_merge.chr{g_count,1}=c_chr;
            gene_merge.start(g_count,1)=min(gene_exons.start(start_g:counter));
            gene_merge.end(g_count,1)=max(gene_exons.end(start_g:counter));
            start_g=counter+1;
            counter=counter+1;
            g_count=g_count+1;
        end
    else
          gene_merge.gene{g_count,1}=c_gene;
            gene_merge.chr{g_count,1}=c_chr;
            gene_merge.start(g_count,1)=min(gene_exons.start(start_g:counter));
            gene_merge.end(g_count,1)=max(gene_exons.end(start_g:counter));
            start_g=counter+1;
            counter=counter+1;
            g_count=g_count+1;
    
    end
end
% CKS1B is not on chr5
% find(ismember(gene_merge.gene,'CKS1B'))
% ans =
%         1290
%        14786
%
hack=ones(slength(gene_merge),1);
hack(14786)=0;
gene_merge=reorder_struct(gene_merge,hack==1);
% 
% find(ismember(gene_merge.gene,'LSP1'))
% ans =
%         2883
%         5177

hack=ones(slength(gene_merge),1);
hack(5177)=0;
gene_merge=reorder_struct(gene_merge,hack==1);

save_struct(gene_merge,'/Users/amaro/Documents/GeneExonMerged_Gencode.tsv')


