%gistic_data 
disp(strcat('Loading Gistic Data'));
CN_Events=load_struct('/Users/amaro/Downloads/HemangiomaPrint/all_thresholded.by_genes_Hemangioma.txt');
%ABSOLUTE maf
disp(strcat('Loading MAF file'));
ABS_maf=load_struct('/Users/amaro/Downloads/HemangiomaPrint/agg_hemangioma_maf.txt');
%sig genes
disp(strcat('Loading Significant Genes'));
sig_genes=load_struct('/Users/amaro/Downloads/cancer_gene_census.tsv');

%ABSOLUTE table
disp(strcat('Loading ABSOLUTE table'));
absolute_table=load_struct('/Users/amaro/Downloads/HemangiomaPrint/Hemangioma_05.07.2014_BlacklistedSegments.ATW.ABSOLUTE.table.txt');





ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.i_judgement,'KEEP'));
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Variant_Classification,'Missense_Mutation')|ismember(ABS_maf.Variant_Classification,'Nonsense_Mutation')|...
    ismember(ABS_maf.Variant_Classification,'Splice_Site')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Ins')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Del'));
samples.names=unique(ABS_maf.sample);
sig_genes=reorder_struct(sig_genes,ismember(sig_genes.gene,ABS_maf.Hugo_Symbol));
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Hugo_Symbol,sig_genes.gene));

NV=slength(sig_genes);
%number of samples
N=length(samples.names);
subplot('position',[0.1 0.05 .77 0.75]) %left bottom width hieght
axis([-1.5 N+2 -2 NV+4]);
set(gca,'visible','off');
%setting dims of boxs for amps and deletions
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;


hom_mut_color=[0 0 0];


for v=1:NV
    for s=1:slength(samples)
        
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end



% place shapes
for v=1:NV
    for s=1:slength(samples)
            %samples.(sig_genes.gene{v})(s,1)=0;

        
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut]; %square
    j=NV-(v-1);
    y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
    
    x2=[s-dx_mut s s+dx_mut s s-dx_mut]; %diamond
    j=NV-(v-1);
    y2=[j j+dy_mut j j-dy_mut j];
    
    x3=[s-dx_mut s+dx_mut s-dx_mut s+dx_mut s-dx_mut]; %double triangle
    y3=[j+dy_mut j+dy_mut j-dy_mut j-dy_mut j+dy_mut];
   
     x4=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
     y4=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    
    
    k=ismember(ABS_maf.sample,samples.names{s});
    
    if ismember(sig_genes.gene{v},ABS_maf.Hugo_Symbol(k))
        c=find(ismember(CN_Events.GeneSymbol,sig_genes.gene{v}));

        % get the position of mutation
        l=find(ismember(ABS_maf.Hugo_Symbol,sig_genes.gene{v})&ismember(ABS_maf.sample,samples.names{s}));
        %samples.(sig_genes.gene{v})(s,1)=1;
        
        %there are mutations check to see if there are scnas too
        
        if ~isempty(samples.gistic_id{s})
        cn=str2double(CN_Events.(samples.gistic_id{s})(c));
        else
            cn=0;
        end
       
        
         if cn<-1
             
             patch(x4,y4,[135/255*.5 206/255*.5 250/255],'edgecolor',[135/255 206/255 250/255])
             patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color)
         elseif cn>1
            patch(x4,y4,[255/255 170/255*.5 135/255*.5],'edgecolor',[255/255 170/255 135/255]);
            patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color,'LineWidth',1.75)
         elseif cn>=-1&&cn<=1
            patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color)
            
            
        end
    else
                c=find(ismember(CN_Events.GeneSymbol,sig_genes.gene{v}));

        % there arent any mutations in this gene
        l=find(ismember(CN_Events.GeneSymbol,sig_genes.gene{v}));
        if ~isempty(samples.gistic_id{s})
        cn=str2double(CN_Events.(samples.gistic_id{s})(c));
        else
            cn=0;
        end
    if cn<-1
        %samples.(sig_genes.gene{v})(s,1)=1;
        patch(x4,y4,[135/255*.5 206/255*.5 250/255],'edgecolor',[135/255 206/255 250/255]);
    end
    if cn>1
        %samples.(sig_genes.gene{v})(s,1)=1;
        patch(x4,y4,[255/255 170/255*.5 135/255*.5],'edgecolor',[255/255 170/255 135/255]);
    end
    
    end
    end
end