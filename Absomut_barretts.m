%% CCF_DP clustering comut 
% need a better name

%ABSOLUTE aggregated maf
disp(strcat('Loading MAF file'));
ABS_maf=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ABSOLUTE.OncoGenes.maf');
%Genes to display
disp(strcat('Loading Significant Genes'));
sig_genes=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/genes_for_ESO.txt');

%ABSOLUTE table
disp(strcat('Loading ABSOLUTE table'));
absolute_table=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/BarrettsABSOLUTEtable.txt');

%Aggregated CCF DP maf
%ccf_dp_maf=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/AggregatedCCFDPmaf.txt');

%sample mapping
%map=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/samplemapping.txt');


disp(strcat('Loading Gene Table'));
Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');

disp(strcat('ABS_seg_file'));

ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Aggregate_Seg_file.txt');

ABS_seg_file.Chromosome=str2double(ABS_seg_file.Chromosome);
ABS_seg_file.Endbp=str2double(ABS_seg_file.Endbp);
ABS_seg_file.Startbp=str2double(ABS_seg_file.Startbp);
ABS_seg_file.total_HZ=str2double(ABS_seg_file.total_HZ);
ABS_seg_file.HZ=str2double(ABS_seg_file.HZ);
ABS_seg_file.SC_HZ=str2double(ABS_seg_file.SC_HZ);
ABS_seg_file.corrected_total_cn=str2double(ABS_seg_file.corrected_total_cn);

samples=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Oncogene_sample_map.txt');
absolute_table=reorder_struct(absolute_table,ismember(absolute_table.array,samples.pair_id));

absolute_table.ploidy=str2double(absolute_table.ploidy);

ABS_maf.ccf=str2double(ABS_maf.ccf_hat);

bottom_alignment=.77;
%% plotting rates
figure('position',[100 100 1000 1200])
hold on



%% main figure
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Hugo_Symbol,sig_genes.gene));
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.i_judgement,'KEEP'));
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Variant_Classification,'Missense_Mutation')|ismember(ABS_maf.Variant_Classification,'Nonsense_Mutation')|...
    ismember(ABS_maf.Variant_Classification,'Splice_Site')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Ins')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Del'));


absolute_table.Genomedoublings=str2double(absolute_table.Genomedoublings);                        
%number of variants
NV=slength(sig_genes);
%number of samples
N=length(samples.names);

%setting axis
%subplot('position',[0.1 0.05 bottom_alignment 0.75]) %left bottom width hieght
axis([-1.5 N+2 -2 NV+4]);
set(gca,'visible','off');
%setting dims of boxs for amps and deletions
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;


hom_mut_color=[0 0 0];


% genome doubling track
j=NV+1;
 for s=1:slength(samples)
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    l=find(ismember(absolute_table.array,samples.pair_id{s}));
    if absolute_table.Genomedoublings(l)==0
            patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])
            samples.GD{s,1}='0';
    else
            patch(x1,y1,[0 0 0],'edgecolor',[0 0 0])
            samples.GD{s,1}='1';
    end
    
   
 end
 
 
  text (-.5,j-.5,'GD','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)

 
 Dys_c=[255/255 182/255 193/255];
 % tissue type track
 j=NV+2;
 for s=1:slength(samples)
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
 
 if isequal(samples.tissue{s},'ND')
  patch(x1,y1,[255/255 182/255 193/255],'edgecolor',[255/255 182/255 193/255])
 elseif isequal(samples.tissue{s},'ESO')
   patch(x1,y1,[34/255 139/255 34/255],'edgecolor',[34/255 139/255 34/255])

 elseif isequal(samples.tissue{s},'Dysplasia')
       patch(x1,y1,[218/255 112/255 214/255],'edgecolor',[218/255 112/255 214/255])

 end
 end
 
   text (-.5,j-.5,'Tissue','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)

 
   
   
 
for v=1:NV
    for s=1:slength(samples)
        
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end
% 

% 



% place shapes
for v=1:NV
    for s=1:slength(samples)
            %samples.(sig_genes.gene{v}){s,1}='0';

        
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
     [abs_cn,HZ]=get_copy_number_of_gene(samples.pair_id{s},ABS_seg_file,sig_genes.gene{v},Gene_Table);
     ploidy=absolute_table.ploidy(ismember(absolute_table.sample,samples.pair_id{s}));
    if ismember(sig_genes.gene{v},ABS_maf.Hugo_Symbol(k))
        
        % get the position of mutation
        l=find(ismember(ABS_maf.Hugo_Symbol,sig_genes.gene{v})&ismember(ABS_maf.sample,samples.names{s}));
        %samples.(sig_genes.gene{v}){s,1}='1';
        
        %there are mutations check to see if there are scnas too
        
      
       if HZ==1 %if cn<-1
        
         
             
             %patch(x4,y4,[135/255*.5 206/255*.5 250/255],'edgecolor',[135/255 206/255 250/255])
             patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color)
       elseif (abs_cn/ploidy)>=2 
            
            if abs_cn>20
            abs_cn=20;
            end
            patch(x4,y4,[255/255 (12*(20-abs_cn))/255 (12*(20-abs_cn))/255],'edgecolor',[255/255 (12*(20-abs_cn))/255 (12*(20-abs_cn))/255]);
                    

           
            patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color,'LineWidth',1.75)
         elseif HZ==0&&(abs_cn/absolute_table.ploidy(ismember(absolute_table.sample,samples.pair_id{s})))<2
            patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color)
            
            
        end
    else
        % there arent any mutations in this gene
        
    if HZ==1
        %samples.(sig_genes.gene{v}){s,1}='1';
      %  patch(x4,y4,[135/255*.5 206/255*.5 250/255],'edgecolor',[135/255 206/255 250/255]);
    end
    if (abs_cn/ploidy)>=2
       %samples.(sig_genes.gene{v}){s,1}='1';
        abs_cn=get_copy_number_of_gene(samples.pair_id{s},ABS_seg_file,sig_genes.gene{v},Gene_Table);
         if abs_cn>20
            abs_cn=20;
            end
            patch(x4,y4,[255/255 (12*(20-abs_cn))/255 (12*(20-abs_cn))/255],'edgecolor',[255/255 (12*(20-abs_cn))/255 (12*(20-abs_cn))/255]);
                    

        
    end
    
    end
    end
end



%labels
for v=1:NV
text (0,NV-(v-.65),sig_genes.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end

for s=1:slength(samples)
h=text(s+dx_mut,0,samples.individual_id{s},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12);
set(h, 'rotation', 90)
end



