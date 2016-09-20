germline_maf='/Users/amaro/Downloads/JMML_LandscapePlot/GermlineLandScapeMafJMML.txt';

% see example file format
G=load_struct(germline_maf);
samples.names=unique(G.Tumor_Sample_Barcode);
% can subsititue for mutsig list 
si.gene={'OUTCOME','NF1','KRAS','NRAS','CBL','PTPN11','RRAS', ...
    'RRAS2', 'SH2B3', 'SETBP1','GATA2', 'RUNX1','ASXL1', 'EZH2', 'DNMT3A','ZRSR2','Mono_7'};
si.gene=si.gene';



N=length(unique(G.Tumor_Sample_Barcode));
NV=length(si.gene);
% code for sorting the plot to show mutual exclusivity
for v=2:NV
    for s=1:N
         k=ismember(G.Tumor_Sample_Barcode,samples.names{s});
         if ismember(si.gene{v},G.Hugo_Symbol(k));
             samples.(si.gene{v})(s,1)=1;
         else
             samples.(si.gene{v})(s,1)=0;
         end
        
    end
end

fields=fieldnames(samples);
samples=sort_struct(samples,{fields{2:end}},(ones(length({fields{2:end}}))*-1));

% default RGB colors
somatic_mut_color=[215/255 48/255 39/255];%[34/255 139/255 34/255];
germ_mut_color=[215/255 48/255 39/255];
miss_germ_mut_color=[69/255 117/255 180/255];
miss_som_mut_color=[69/255 117/255 180/255]; %[116/255 173/255 209/255];
deletion_color=[69/255 117/255 180/255];

hold on


%N=100%length(unique(G.sample));
%NV=20%slength(si);
axis([-1.5 N+2 -2 NV+8]);
set(gca,'visible','off');
%setting dims of boxs 
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;

% placing all the blank panels
for v=1:NV
    for s=1:N
     %these correspond to the shapes of the panels in this case rectangles   
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1)+4;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.95 .95 .95],'edgecolor',[.95 .95 .95])

        
    end
end

% place germline / somatic / LOH symbols
for v=2:NV
    for s=1:N
        % loop through the variants and samples placing patches for
        % mutations in either the germline half or somatic half based on
        % missense and nonsense
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];%
     j=NV-(v-1)+4;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    
    y_germ=[j-(dy_cn-.005) j-(dy_cn-.005) j j j-dy_cn];% germline y
    y_som=[j j j+(dy_cn-.005) j+(dy_cn-.005) j]; %somatic y
    
    k=ismember(G.Tumor_Sample_Barcode,samples.names{s});
    if ismember(si.gene{v},G.Hugo_Symbol(k));
        l=find(ismember(G.Hugo_Symbol,si.gene{v})&ismember(G.Tumor_Sample_Barcode,samples.names{s}));
        for i=1:length(l)
           
            if isequal(G.is_germline{l(i)},'0')&&~(isequal(G.is_LOH{l(i)},'2'))
                if isequal(G.Variant_Class{l(i)},'Nonsense')
                    if isequal(G.Relapse{l(i)},'0')
                patch(x1,y_som,somatic_mut_color,'edgecolor',somatic_mut_color,'EdgeAlpha',0)
                    else
                patch(x1,y_som,somatic_mut_color,'edgecolor',somatic_mut_color,'EdgeAlpha',0)
                    for ff=1:2:8
                        x2=[(s-dx_mut)+.1*(ff-1) (s-dx_mut)+.1*(ff) (s-dx_mut)+.1*(ff) (s-dx_mut)+.1*(ff-1) (s-dx_mut)+.1*(ff-1)];
                        patch(x2,y_som,[0 0 0],'edgecolor',[0 0 0],'EdgeAlpha',0,'LineWidth',.25)
                    end
                    end
                elseif isequal(G.Relapse{l(i)},'0')
                    patch(x1,y_som,miss_som_mut_color,'edgecolor',miss_som_mut_color,'EdgeAlpha',0)
                else
                    
                    patch(x1,y_som,miss_som_mut_color,'edgecolor',[0 0 0],'EdgeAlpha',0)
                    
                    for ff=1:2:8
                        x2=[(s-dx_mut)+.1*(ff-1) (s-dx_mut)+.1*(ff) (s-dx_mut)+.1*(ff) (s-dx_mut)+.1*(ff-1) (s-dx_mut)+.1*(ff-1)];
                        patch(x2,y_som,[0 0 0],'edgecolor',[0 0 0],'EdgeAlpha',0,'LineWidth',.25)
                    end
                end
            elseif ~(isequal(G.is_LOH{l(i)},'2'))
                if isequal(G.Variant_Class{l(i)},'Nonsense')
                patch(x1,y_germ,germ_mut_color,'edgecolor',germ_mut_color,'EdgeAlpha',0)
                else
                    patch(x1,y_germ,miss_germ_mut_color,'edgecolor',miss_germ_mut_color,'EdgeAlpha',0)
                end
            end
            
             if isequal(G.is_LOH{l(i)},'1')
                patch(x1,y1,[.95 .95 .95],'edgecolor',[77/255 77/255 77/255],'LineWidth',1,'FaceAlpha',0)

            elseif isequal(G.is_LOH{l(i)},'2')
                
                patch(x1,y_som,[69/255 117/255 180/255],'edgecolor',[69/255 117/255 180/255],'LineWidth',2)
               
            end
            
        
        end
            
        
    end
        

    end
end

v=1;
for s=1:N
        
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];%
    j=NV-(v-1)+4;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    k=ismember(G.Tumor_Sample_Barcode,samples.names{s});
    if ismember(si.gene{v},G.Hugo_Symbol(k));
        l=find(ismember(G.Hugo_Symbol,si.gene{v})&ismember(G.Tumor_Sample_Barcode,samples.names{s}));
    if isequal(G.Outcomes{l(i)},'0')
        patch(x1,y1,[.8 .8 .8],'edgecolor',[.8 .8 .8],'LineWidth',2)
    elseif isequal(G.Outcomes{l(i)},'1')
        patch(x1,y1,[.2 .2 .2],'edgecolor',[.2 .2 .2],'LineWidth',2)
    elseif isequal(G.Outcomes{l(i)},'3')
        patch(x1,y1,[1 1 1],'edgecolor',[1 1 1],'LineWidth',2)
    end
    end
end
% extra track for monosomy 7 this is an easy place to put more tracks. 
si.gene{18,1}='Mono-7';

si.gene{1,1}='Outcome';


for v=1:NV
    if v>16
        text (.5,NV-(v-1)+3.5,si.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)

    else
text (.5,NV-(v-1)+3.5,si.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
    end
end


info=reorder_struct(G,ismember(G.Hugo_Symbol,'INFO'));
si.info={'Gender';'WBC';'Platelet';'Age'};

for v=1:length(si.info)
    j=NV-(NV-length(si.gene)-(v))+3;
    text (.5,j-.5,si.info{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end

for v=1:length(si.info)
    for s=1:N
        % loop through the variants and samples placing patches for
        % mutations in either the germline half or somatic half based on
        % missense and nonsense
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];%
    j=NV-(NV-length(si.gene)-(v))+3;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    k=find(ismember(info.Tumor_Sample_Barcode,samples.names{s}));
    if isequal(si.info{v},'Gender')
        if isequal(info.Gender{k},'1')
            patch(x1,y1,[100/255 100/255 255/255],'edgecolor',[100/255 100/255 255/255],'LineWidth',2)
        elseif isequal(info.Gender{k},'2')
            patch(x1,y1,[1 105/255 180/255],'edgecolor',[1 105/255 180/255],'LineWidth',2)
        end
    end
    
    if isequal(si.info{v},'Age')
        age=str2double(info.Age{k})/96;
        if age>0
        y1=[j-dy_cn j-dy_cn j-dy_cn+(2*dy_cn*age) j+-dy_cn+(2*dy_cn*age) j-dy_cn];
            patch(x1,y1,[0 0 0],'edgecolor',[0 0 0],'LineWidth',2)   
        end
    end
    
    if isequal(si.info{v},'WBC')
        age=str2double(info.WBC{k})/200;
        if age>0
        y1=[j-dy_cn j-dy_cn j-dy_cn+(2*dy_cn*age) j+-dy_cn+(2*dy_cn*age) j-dy_cn];
            patch(x1,y1,[0 0 0],'edgecolor',[0 0 0],'LineWidth',2)   
        end
    end
    
    if isequal(si.info{v},'Platelet')
        age=str2double(info.Platelet{k})/220;
        if age>0
        y1=[j-dy_cn j-dy_cn j-dy_cn+(2*dy_cn*age) j+-dy_cn+(2*dy_cn*age) j-dy_cn];
            patch(x1,y1,[0 0 0],'edgecolor',[0 0 0],'LineWidth',2)   
        end
    end
    
    end
end

% legends!


% j=-1;
% s=1;
% x_1=[s-dx_mut s s s-dx_mut s-dx_mut];
% x_2=[s s+dx_mut s+dx_mut s s];
% y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
% text(1+.5,-1.5,'Germline Mut.','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
% patch(x_1,y1,germ_mut_color,'edgecolor',germ_mut_color)
% patch(x_2,y1,miss_germ_mut_color,'edgecolor',miss_germ_mut_color)
% 
% j=-1;
% s=9;
% x_1=[s-dx_mut s s s-dx_mut s-dx_mut];
% x_2=[s s+dx_mut s+dx_mut s s];
% y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
% text(s+.5,-1.5,'Somatic Mut.','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
% patch(x_1,y1,somatic_mut_color,'edgecolor',somatic_mut_color)
% patch(x_2,y1,miss_som_mut_color,'edgecolor',miss_som_mut_color)
% 
% 
% s=17;
% x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
% y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
% text(s+.5,-1.5,'Biallelic.','VerticalAlignment','bottom','HorizontalAlignment','Left','FontSize',16);
% patch(x1,y1,[.9 .9 .9],'edgecolor',deletion_color)
% 
% 
% 

for s=1:N
h=text(s+dx_mut,4,samples.names{s},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12);
set(h, 'rotation', 90)
end
