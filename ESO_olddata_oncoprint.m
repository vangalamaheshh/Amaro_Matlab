%ESO old data onco print

%gistic_data 
disp(strcat('Loading Gistic Data'));
%CN_Events=load_struct('/Users/amaro/Downloads/OLD_Eso_oncoprint/Old_ESO_all_thresholded.by_genes.txt');
%ABSOLUTE maf
disp(strcat('Loading MAF file'));
ABS_maf=load_struct('/Users/amaro/Downloads/OLD_Eso_Oncoprint/ESOmaffilt.txt');%sig genes
disp(strcat('Loading Significant Genes'));
sig_genes=load_struct('/Users/amaro/Downloads/BarrettsCCFDp/genesForESO_144_oncogenes.txt');

%ABSOLUTE table
disp(strcat('Loading ABSOLUTE table'));
absolute_table=load_struct('/Users/amaro/Downloads/OLD_Eso_Oncoprint/146_samples.ATW.ABSOLUTE.table.txt');
disp(strcat('Loading Gene Table'));
Gene_Table=load_struct('/Users/amaro/Downloads/RefSeqGenes.txt');
ABS_seg_file=load_struct('/Users/amaro/Downloads/OLD_Eso_Oncoprint/Agg_SegFile.txt');



ABS_seg_file.Chromosome=str2double(ABS_seg_file.Chromosome);
ABS_seg_file.Endbp=str2double(ABS_seg_file.Endbp);
ABS_seg_file.Startbp=str2double(ABS_seg_file.Startbp);
ABS_seg_file.total_HZ=str2double(ABS_seg_file.total_HZ);
ABS_seg_file.HZ=str2double(ABS_seg_file.HZ);
ABS_seg_file.SC_HZ=str2double(ABS_seg_file.SC_HZ);
ABS_seg_file.corrected_total_cn=str2double(ABS_seg_file.corrected_total_cn);

% ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Hugo_Symbol,sig_genes.gene));
% ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.i_judgement,'KEEP'));
% ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Variant_Classification,'Missense_Mutation')|ismember(ABS_maf.Variant_Classification,'Nonsense_Mutation')|...
%     ismember(ABS_maf.Variant_Classification,'Splice_Site')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Ins')...
%     |ismember(ABS_maf.Variant_Classification,'Frame_Shift_Del')|ismember(ABS_maf.Variant_Classification,'In_Frame_Del')...
%     |ismember(ABS_maf.Variant_Classification,'In_Frame_Ins'));

%subset events to sig_genes
%CN_Events=reorder_struct(CN_Events,ismember(CN_Events.GeneSymbol,sig_genes.gene));
% determine mutation rate    %ESO specific sample name fix
%                         for i=1:slength(ABS_maf)
%                             strs=regexp(ABS_maf.sample{i},'-','split');
%                             ABS_maf.sample{i}=strcat(strs{1},strs{2},'TPN');
%                         end
%                         for i=1:slength(ABS_seg_file)
%                             strs=regexp(ABS_seg_file.sample{i},'-','split');
%                             ABS_seg_file.sample{i}=strcat(strs{1},strs{2},'TPN');
%                         end
%                         

                        
 samples_ord=load_struct('ESO_Samples_Events.txt'); 
 samples.names=samples_ord.names;
 samples.pair_id=samples_ord.pair_id;
absolute_table.Genomedoublings=str2double(absolute_table.Genomedoublings);                        

                        %ESO specific sample name fix
%                         for i=1:slength(absolute_table)
%                          strs=regexp(absolute_table.array{i},'-','split');
%                          absolute_table.array{i}=strcat(strs{1},strs{2},'TPN');
%                         end
                        
%number of variants in print equal to number of genes
NV=slength(sig_genes);
%number of samples:
%CN_fields=fieldnames(CN_Events);
%N=length(CN_fields)-3;
N=slength(samples_ord);
%%turning copynumbers to doubles
%for i=3:length(CN_fields)
%CN_Events.(CN_fields{i})=str2double(CN_Events.(CN_fields{i}));
%end


%sorting code
samples.GD=str2double(samples_ord.GD);

for i=1:slength(sig_genes)
    
   
    for j=1:slength(samples)
%        dl=find(ismember(absolute_table.array,samples.names{j}));
%       % l=find(ismember(CN_Events.GeneSymbol,sig_genes.gene{i}));
%        k=ismember(ABS_maf.sample,samples.names{j});
%     if absolute_table.Genomedoublings(dl)>0
%         samples.doubled(j,1)=1;
%     else
%         samples.doubled(j,1)=0;
%     end
        
        
        
        if isequal(samples_ord.(sig_genes.gene{i}){j,1},'None')
        samples.(sig_genes.gene{i})(j,1)=0;
        else
            samples.(sig_genes.gene{i})(j,1)=1;
        end
    
    
        
        
       
    end
    
end
samples.GD(samples.GD>1)=1;
f_samples=fieldnames(samples);
samples=sort_struct(samples,{f_samples{3:end}},ones(length(f_samples)-2,1)*-1);


% %% rate barplot
% for s=1:slength(samples)
% rate.numberofClonalmutations(s,1)=sum(ismember(ABS_maf.sample,samples.names{s})&ismember(ABS_maf.clonalix,'TRUE'));
% rate.numberofSubClonalmutations(s,1)=sum(ismember(ABS_maf.sample,samples.names{s})&ismember(ABS_maf.clonalix,'FALSE'));
% end
% rate.mMbClonal=((rate.numberofClonalmutations)./30000000)*1000000;
% rate.mMbSubClonal=((rate.numberofSubClonalmutations)./30000000)*1000000;
% 

% %plotting rates
% figure('position',[100 100 1000 1200])
% hold on
% 
% subplot('position',[0.125 0.8 0.78 0.15]) %left bottom width hieght
% P=bar(1:slength(samples), [rate.mMbClonal rate.mMbSubClonal],'stack');
% set(P(1),'facecolor',[255/255 127/255 80/255],'EdgeColor',[1,1,1])
% set(P(2),'facecolor',[80/255 127/255 255/255],'EdgeColor',[1,1,1])
% set(gca,'xcolor',[1 1 1],'xtick',[])
% ylabel('Mut per mB'); box off
% limy=get(gca,'ylim');
% limx=get(gca,'xlim');
% Pbaseline=get(P,'BaseLine');
% set(Pbaseline{1}(1),'Color',[.9 .9 .9],'LineWidth',1);
% 
% 
% text(limx(2)-3,limy(2)-1,'Clonal','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',12);
% text(limx(2)-3,limy(2)+2,'Subclonal','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',12);
% j=limy(2)+1; s=limx(2)-4;
% x1=[s-.75 s+.75 s+.75 s-.75 s-.75];
% y1=[j-.75 j-.75 j+.75 j+.75 j-.75];
% patch(x1,y1,[255/255 127/255 80/255],'edgecolor',[255/255 127/255 80/255])
% 
% j=limy(2)+4; s=limx(2)-4;
% x1=[s-.75 s+.75 s+.75 s-.75 s-.75];
% y1=[j-.75 j-.75 j+.75 j+.75 j-.75];
% patch(x1,y1,[80/255 127/255 255/255],'edgecolor',[80/255 127/255 255/255])


%% main plot




%setting axis
subplot('position',[0.1 0.05 0.85 0.7]) %left bottom width hieght
axis([-5 N+3 -2 NV+3]);
set(gca,'visible','off');
%setting dims of boxs for amps and deletions
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;







disp(strcat('Plotting'));


deletion_color=[5,112,176]./255;
amp_color=[255/255 170/255 135/255];
het_mut_color=[34/255 139/255 34/255];
hom_mut_color=[34/255 139/255 34/255];


j=-1;
s=1;
x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
text(2+1+dx_mut,-1.15,'Hom. Mut','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
patch(x1,y1,hom_mut_color,'edgecolor',hom_mut_color)

s=21;
x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
text(s+1+dx_mut,-1.15,'Het. Mut','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
patch(x1,y1,het_mut_color,'edgecolor',het_mut_color)


s=41;
x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn]; 
text(s+1+dx_mut,-1.15,'Deletion','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
patch(x1,y1,deletion_color,'edgecolor',[.9 .9 .9]);


s=61;
x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn]; 
text(s+1+dx_mut,-1.15,'Amplification','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
patch(x1,y1,[1 abs(1/2) abs(1/2)].*amp_color,'edgecolor',[.9 .9 .9]);

% genome doublings track
j=NV+2;
 for s=1:slength(samples)
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
     l=find(ismember(absolute_table.array,samples.names{s}));
    if absolute_table.Genomedoublings(l)==0
           patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])
    else
            patch(x1,y1,[0 0 0],'edgecolor',[0 0 0])
    end
    
   
 end

 text (-.5,NV+1.3,'GD','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
 
 
 
%place all the patches

for v=1:NV
    for s=1:slength(samples)
        
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end


ss=zeros(NV,length(samples));

%place the mutations



%place all the copynumber events
absolute_table.ploidy=str2double(absolute_table.ploidy);
for v=1:NV
    for s=1:slength(samples)
    
        
        
    [abs_cn,HZ]=get_copy_number_of_gene(samples.names{s},ABS_seg_file,sig_genes.gene{v},Gene_Table);
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];   
    
    %l=find(ismember(CN_Events.GeneSymbol,sig_genes.gene{v}));
    %cn=CN_Events.(samples.names{s})(l);
%     if HZ==1
%         patch(x1,y1,deletion_color,'edgecolor',deletion_color);
%         ss(v,s)=1;
%     end
   if (abs_cn/absolute_table.ploidy(ismember(absolute_table.array,samples.names{s})))>=2
        
        
            if abs_cn>20
            abs_cn=20;
            end
            patch(x1,y1,[255/255 (12*(20-abs_cn))/255 (12*(20-abs_cn))/255],'edgecolor',[255/255 (12*(20-abs_cn))/255 (12*(20-abs_cn))/255]);
%            patch(x4,y4,[255/255 (36*(8-abs_cn))/255 (36*(8-abs_cn))/255],'edgecolor',[255/255 (36*(8-abs_cn))/255 (36*(8-abs_cn))/255]);
    end
    
    end
end


% missense=[51,160,44]./255; %green
% nonsense= [227,26,28]./255; %red
% splice= [255,127,0]./255; %orange
% indel= [106,61,154]./255; %purple
missense=[0 0 0];
nonsense=[0 0 0];
splice=[0 0 0];
indel=[0 0 0];
for v=1:NV
    for s=1:slength(samples)
    
        
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
    j=NV-(v-1);
    y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
    
    k=ismember(ABS_maf.sample,samples.names{s});
    
    if ismember(sig_genes.gene{v},ABS_maf.Hugo_Symbol(k))
     l=find(ismember(ABS_maf.Hugo_Symbol,sig_genes.gene{v})&ismember(ABS_maf.sample,samples.names{s}));
     l=l(1);
   x2=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y2=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    %patch(x2,y2,[.9 .9 .9],'edgecolor',[.9 .9 .9])
    if isequal(ABS_maf.Variant_Type{l},'DEL')||isequal(ABS_maf.Variant_Type{l},'INS')
          patch(x1,y1,indel,'edgecolor',indel)
    elseif isequal(ABS_maf.Variant_Classification{l},'Missense_Mutation')
                  patch(x1,y1,missense,'edgecolor',missense)
    elseif isequal(ABS_maf.Variant_Classification{l},'Nonsense_Mutation')
                  patch(x1,y1,nonsense,'edgecolor',nonsense)
    elseif isequal(ABS_maf.Variant_Classification{l},'Splice_Site')
                  patch(x1,y1,splice,'edgecolor',splice)
    end
                    ss(v,s)=1;

        %end
    end
    end
end



%labels
for v=1:NV
text (-.5,NV-(v-.7),sig_genes.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end




