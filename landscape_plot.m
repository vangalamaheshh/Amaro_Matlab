function landscape_plot(varargin)
% Firehose Landscape Plot Module
% get input files
files=options_parsing_for_landscape_plot(varargin);


if isempty(files.sample_to_pair)
    error('Sample to Pair mapping is required for plotting')
end

if isempty(files.ABS_maf)&&isempty(files.Somatic_maf)&&isempty(files.Germ_maf)
    error('At least one maf is required')
end

if isempty(files.sig_genes)
    error('list of genes to display required')
end

if ~isempty(files.gistic_table)
 
 files.gistic_table=reorder_struct(files.gistic_table,ismember(files.gistic_table.GeneSymbol,files.sig_genes.gene));
 
 CN_fields=fieldnames(files.gistic_table);
 for i=3:length(CN_fields)
    files.gistic_table.(CN_fields{i})=str2double(files.gistic_table.(CN_fields{i}));
 end
 for i=1:slength(files.sample_to_pair)
  strs=regexp(files.sample_to_pair.pair_id{i},'-','split');
  files.sample_to_pair.gistic_id{i,1}=strcat(strs{:});
  % This is a hack for old capseg naming bugs to be removed as soon as
  % recapseg becomes common
  if isequal(files.sample_to_pair.gistic_id{i,1}(end),'T')
    files.sample_to_pair.gistic_id{i,1}=files.sample_to_pair.gistic_id{i,1}(1:end-1);
  end
 end
 
end


NV=slength(files.sig_genes);
N=slength(files.sample_to_pair);


%subplot('position',[0.1 0.05 .77 0.75]) %left bottom width hieght
axis([-1.5 N+2 -2 NV+4]);
set(gca,'visible','off');
%setting dims of boxs for amps and deletions
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;
dx_g=.1; dy_g=.1;
mut_c=[0 100/255 0];
G_c=[255/255 215/255 32/255];
deletion_c=[0/255 0/255 225/255];
amplification_c=[255/255 0/255 0/255];



% if absolute table is passed fill in a track for Genome Doublings
 if ~isempty(files.absolute_table)
     disp('Plotting Genome Doublings')
    Genomedoublings=str2double(files.absolute_table.Genomedoublings);                        
    j=NV+1;
    for s=1:slength(files.sample_to_pair)
        x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
        l=find(ismember(files.absolute_table.array,files.sample_to_pair.pair_id{s}));
    if Genomedoublings(l)==0
            patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])
    else
            patch(x1,y1,[0 0 0],'edgecolor',[0 0 0])
    end
    end
    text (-.5,j-.5,'GD','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
 end


 % place base line patches
 for v=1:NV
     
    for s=1:slength(files.sample_to_pair)
        
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])  
    end
 end
 
 
% if gistic table is passed place CN patches
if ~isempty(files.gistic_table)
    disp('Plotting CN panels')
    for v=1:NV
    for s=1:slength(files.sample_to_pair)
    
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];   
    
    l=find(ismember(files.gistic_table.GeneSymbol,files.sig_genes.gene{v}));
    cn=files.gistic_table.(files.sample_to_pair.gistic_id{s})(l);
    if cn<-1
        patch(x1,y1,deletion_c,'edgecolor',[.9 .9 .9]);
    end
    if cn>1
        patch(x1,y1,amplification_c,'edgecolor',[.9 .9 .9]);
    end
    
    end
    end
end


% if there is an absolute maf perfer that one to the somatic one
% if there is a germline maf we combine that with the somatic / absolute
% one to display double hits. 

if ~isempty(files.Germ_maf) && ~isempty(files.ABS_maf)
    s_maf='ABS_maf';
    g_maf=1;
    
elseif ~isempty(files.Germ_maf)&&~isempty(files.Somatic_maf)
    s_maf='Somatic_maf';
    g_maf=1;
    
elseif ~isempty(files.ABS_maf)
    s_maf='ABS_maf';
    g_maf=0;

elseif ~isempty(files.Somatic_maf)
    s_maf='Somatic_maf';
    g_maf=0;
   

end

files.(s_maf)=reorder_struct(files.(s_maf),ismember(files.(s_maf).Variant_Classification,'Missense_Mutation')|ismember(files.(s_maf).Variant_Classification,'Nonsense_Mutation')|...
    ismember(files.(s_maf).Variant_Classification,'Splice_Site')|ismember(files.(s_maf).Variant_Classification,'Frame_Shift_Ins')|ismember(files.(s_maf).Variant_Classification,'Frame_Shift_Del'));
files.(s_maf)=reorder_struct(files.(s_maf),ismember(files.(s_maf).i_judgement,'KEEP'));

for v=1:NV
    for s=1:slength(files.sample_to_pair)   
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
    j=NV-(v-1);
    y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
    
    k=ismember(files.(s_maf).Tumor_Sample_UUID,files.sample_to_pair.sample_id{s});
    
    if ismember(files.sig_genes.gene{v},files.(s_maf).Hugo_Symbol(k))
        
    l=find(ismember(files.(s_maf).Hugo_Symbol,files.sig_genes.gene{v})&ismember(files.(s_maf).Tumor_Sample_UUID,files.sample_to_pair.sample_id{s}));
    patch(x1,y1,mut_c,'edgecolor',mut_c)    
    end        
    end
end

if g_maf
    for v=1:NV
        for s=1:slength(files.sample_to_pair)
          x1=[s-dx_g s+dx_g s+dx_g s-dx_g s-dx_g];
          j=NV-(v-1);
          y1=[j-dy_g j-dy_g j+dy_g j+dy_g j-dy_g]; 
          k=ismember(files.Germ_maf.Tumor_Sample_UUID,files.sample_to_pair.sample_id{s});
          if ismember(files.sig_genes.gene{v},files.Germ_maf.Hugo_Symbol(k))
            patch(x1,y1,G_c,'edgecolor',G_c)  
          end
        
        end
    end

end
%labels
for v=1:NV
text (0,NV-(v-.65),files.sig_genes.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end





end

