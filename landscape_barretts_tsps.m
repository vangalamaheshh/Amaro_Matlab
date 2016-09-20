%case 27 and 109 CDKN2A not shared
% related  shifted up
% red box around stray sample. 
% 63 needs red box

%ABSOLUTE aggregated maf
disp(strcat('Loading MAF file'));
ABS_maf=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/TSPsMaf.txt');
%Genes to display
disp(strcat('Loading Significant Genes'));
sig_genes=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/genes_for_ESO_TSPs.txt');

%ABSOLUTE table
disp(strcat('Loading ABSOLUTE table'));
absolute_table=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/BarrettsABSOLUTEtable.txt');

%Aggregated CCF DP maf
%ccf_dp_maf=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/AggregatedCCFDPmaf.txt');

%sample mapping
map=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/samplemapping.txt');
old_sample_table=load_struct('Barretts_Samples_Events_onc_new_order.txt');
map=reorder_struct(map,ismember(map.pair_id,old_sample_table.pair_id));
%gistic_data 
disp(strcat('Loading Gistic Data'));
CN_Events=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/all_thresholded.by_genes_Barretts.txt');

CN_Events=reorder_struct(CN_Events,ismember(CN_Events.GeneSymbol,sig_genes.gene));

ABS_maf.position=ABS_maf.Start_position;

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

CN_fields=fieldnames(CN_Events);
for i=3:length(CN_fields)
CN_Events.(CN_fields{i})=str2double(CN_Events.(CN_fields{i}));
end


samples=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/TSPs_Sample_Map.txt');
%samples=sort_struct(samples,'relatedness');
samples.Individual=samples.individual_id;
samples.names=samples.case_sample;
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.sample,samples.case_sample));

ABS_maf.ccf=str2double(ABS_maf.ccf_hat);
figure('position',[100 100 1000 1200])
hold on
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Hugo_Symbol,sig_genes.gene));
ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.i_judgement,'KEEP'));

absolute_table.Genomedoublings=str2double(absolute_table.Genomedoublings);                        
%number of variants
NV=slength(sig_genes);
%number of samples
N=length(unique(samples.Individual));
load('Is.mat');
%setting axis
%subplot('position',[0.1 0.05 bottom_alignment 0.75]) %left bottom width hieght
axis([-1.5 N+2 -2 NV+4]);
set(gca,'visible','off');
%setting dims of boxs for amps and deletions
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;

b_clo=[218/255 112/255 214/255];
shared=[242/255 125/255 7/255];
e_clo=[34/255 139/255 34/255];
conver=[255/255 0/255 0/255];

for i=1:slength(ABS_maf)
k=ismember(samples.names,ABS_maf.sample{i});
if ~isempty(find(k==1))
ABS_maf.tissue{i,1}=samples.tissue{k};
ABS_maf.Ind{i,1}=samples.Individual{k};
else
    ABS_maf.tissue{i,1}='NA';
    ABS_maf.Ind{i,1}='NA';
end
end


for v=1:NV
    for s=1:N
        
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end


for v=1:NV
    for s=1:N
        
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];%
    
    j=NV-(v-1);
    
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    
    y_b=[j-(dy_cn-.005) j-(dy_cn-.005) j-(dy_cn-.005) j+(dy_cn-.005) j-(dy_cn-.005)];% Barretts
    y_e=[j+(dy_cn-.005) j-(dy_cn-.005) j+(dy_cn-.005) j+(dy_cn-.005) j+(dy_cn-.005)]; %ESO 
    
    k=ismember(ABS_maf.Ind,Is{s});
    
    if ismember(sig_genes.gene{v},ABS_maf.Hugo_Symbol(k));
        l=find(ismember(ABS_maf.Hugo_Symbol,sig_genes.gene{v})&ismember(ABS_maf.Ind,Is{s}));
        for i=1:length(l)
            if isequal(ABS_maf.tissue{l(i)},'Dysplasia')
                
                patch(x1,y_b,b_clo,'edgecolor',b_clo,'EdgeAlpha',0)
            elseif isequal(ABS_maf.tissue{l(i)},'ND')
                    patch(x1,y_b,b_clo,'edgecolor',b_clo,'EdgeAlpha',0)
                
           
            else
                if isequal(ABS_maf.clonalix{l(i)},'TRUE')
                patch(x1,y_e,e_clo,'edgecolor',e_clo,'EdgeAlpha',0)
                else
                    patch(x1,y_e,e_clo,'edgecolor',e_clo,'EdgeAlpha',0)
                end
            end
       
         
            
        end
            
        if length(l)>1
            if ~isequal(ABS_maf.position{l(1)},ABS_maf.position{l(2)})||~(isequal(ABS_maf.Tumor_Seq_Allele2{l(1)},ABS_maf.Tumor_Seq_Allele2{l(2)}))
            %patch(x1,y1,[1 1 1],'edgecolor',conver,'FaceAlpha',0,'LineWidth',1.5)
            else
                patch(x1,y1,shared,'edgecolor',shared,'LineWidth',1.5)
            end
        end
    end
    
    %g=samples.gistic_id(ismember(samples.Individual,Is{s}));
    p=samples.pair_id(ismember(samples.Individual,Is{s}));
    for i=1:length(p)
        %c=find(ismember(CN_Events.GeneSymbol,sig_genes.gene{v}));
        %cn=CN_Events.(p{i})(c);
        [abs_cn,HZ]=get_copy_number_of_gene(p{i},ABS_seg_file,sig_genes.gene{v},Gene_Table);
        if isequal(samples.tissue(ismember(samples.pair_id,p{i})),{'ESO'}) && HZ==1
            patch(x1,y_e,e_clo,'edgecolor',e_clo,'EdgeAlpha',0)
        elseif HZ==1 && ~isequal(samples.tissue(ismember(samples.pair_id,p{i})),{'ESO'})
             patch(x1,y_b,b_clo,'edgecolor',b_clo,'EdgeAlpha',0)
        end
        
%         if isequal(sig_genes.gene{v},'CDKN2A')&& isequal(Is{s},'UMBEER_63')
%              patch(x1,y1,[1 1 1],'edgecolor',conver,'FaceAlpha',0,'LineWidth',1.5)
% 
%         end
%         
            
    
        if isequal(sig_genes.gene{v},'CDKN2A')&& isequal(Is{s},'UMBEER_49')
                patch(x1,y1,shared,'edgecolor',shared,'LineWidth',1.5)

        end
    if isequal(sig_genes.gene{v},'ARID1A')&& isequal(Is{s},'UMBEER_71')
                patch(x1,y1,shared,'edgecolor',shared,'LineWidth',1.5)
    end
                if isequal(sig_genes.gene{v},'TP53')&& isequal(Is{s},'UMBEER_33')
                patch(x1,y1,shared,'edgecolor',shared,'LineWidth',1.5)
                 end
        
        
    end
    end
end

 samples.relatedness=str2double(samples.relatedness);
% relatedness track
 j=NV+3;
 for s=1:N
     
    R=mean(samples.relatedness(ismember(samples.Individual,Is{s})));
   
     x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
     y1=[j-dy_cn j-dy_cn j-dy_cn+(2*dy_cn)*((R/100)) j-dy_cn+(2*dy_cn)*((R/100)) j-dy_cn];

    if R>0
     
    patch(x1,y1,[100/255 149/255 237/255]...
        ,'edgecolor',[100/255 149/255 237/255])
    else
    end
 end
 text (-.5,j-.4,'% Relate','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
 
  j=NV+2;
  s=1;
  s1=1;
 while s1<=N
 x1=[s1-dx_mut s1+dx_mut s1+dx_mut s1-dx_mut s1-dx_mut];%

 y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
 if isequal(samples.tissue{s},'Dysplasia')
     patch(x1,y1,b_clo,'edgecolor',b_clo,'EdgeAlpha',0)
     s=s+1;
     s1=s1+1;
elseif isequal(samples.tissue{s},'ND')
    patch(x1,y1,[255/255 182/255 193/255],'edgecolor',[255/255 182/255 193/255],'EdgeAlpha',0)
    s=s+1;
    s1=s1+1;
 else
     s=s+1;
 end

end
 
 
 
 j=NV+1;
  for s=1:N
     
    gd_cancer=samples.GD(ismember(samples.Individual,Is{s})&ismember(samples.tissue,'ESO'));
    gd_barretts=samples.GD(ismember(samples.Individual,Is{s})&~ismember(samples.tissue,'ESO'));
     x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];%
      y_b=[j-(dy_cn-.005) j-(dy_cn-.005) j-(dy_cn-.005) j+(dy_cn-.005) j-(dy_cn-.005)];% Barretts
    y_e=[j+(dy_cn-.005) j-(dy_cn-.005) j+(dy_cn-.005) j+(dy_cn-.005) j+(dy_cn-.005)]; %ESO 
    patch(x1,y_b,[.9 .9 .9],'edgecolor',[.9 .9 .9])
       patch(x1,y_e,[.9 .9 .9],'edgecolor',[.9 .9 .9])
    if isequal(gd_cancer,{'1'})
    patch(x1,y_e,e_clo,'edgecolor',e_clo,'EdgeAlpha',0)
    end
    if isequal(gd_barretts,{'1'})&&isequal(samples.tissue(ismember(samples.Individual,Is{s})&~ismember(samples.tissue,'ESO')),{'ND'})
        patch(x1,y_b,[255/255 182/255 193/255],'edgecolor',[255/255 182/255 193/255],'EdgeAlpha',0)
    end
    if isequal(gd_barretts,{'1'})&&isequal(samples.tissue(ismember(samples.Individual,Is{s})&~ismember(samples.tissue,'ESO')),{'Dysplasia'})
         patch(x1,y_b,b_clo,'edgecolor',b_clo,'EdgeAlpha',0)
    end
       

  end
 
 
   text (.5,j-.4,'GD','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
   
%    j=NV+2;

%  for s=1:slength(samples)
%     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
%     y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
%     l=samples.pair_id(
%     R=mean(samples.(ismember(samples.Individual,Is{s})));
%     if absolute_table.Genomedoublings(l)==0
%             patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])
%             samples.GD{s,1}='0';
%     else
%             patch(x1,y1,[0 0 0],'edgecolor',[0 0 0])
%             samples.GD{s,1}='1';
%     end
%     
%    
%  end
%   text (-.5,j-.5,'GD','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)

   
   
   


for v=1:NV
text (.5,NV-(v-.8),sig_genes.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end


for s=1:N
h=text(s+dx_mut,0,Is{s},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12);
set(h, 'rotation', 90)
end





