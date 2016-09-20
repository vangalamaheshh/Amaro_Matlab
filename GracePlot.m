addpath('/cga/tcga-gsc/home/lawrence/cga/trunk/matlab/mike')
addpath('/cga/tcga-gsc/home/lawrence/cga/trunk/matlab/')

%M=load_struct('/Users/amaro/Documents/JMML_AllelicCapseg_Plots/JMMLHitsByPatient2.txt');
M=load_struct('~/Documents/somatic_germline_genes.txt');

M.Patient=M.Tumor_Sample_Barcode
M.Gene=M.Hugo_Symbol
M.JustInRelapse=M.Germline

%genes={'NF1','KRAS','NRAS','CBL','PTPN11','RRAS', ...
%    'RRAS2', 'SH2B3','JAK3', 'SETBP1','GATA2', 'RUNX1','ASXL1', 'EZH2', 'DNMT3A','ZRSR2','Mono_7'};
%genes={'ATM', 'BCOR', 'CHD2', 'DDX3X', 'EGR2', 'FBXW7', 'HIST1H1E', 'ITPKB', 'KRAS', 'MED12', 'MYD88', 'NOTCH1', 'NRAS', 'POT1', 'RIPK1', 'SAMHD1', 'SF3B1', 'TP53', 'XPO1', 'ZMYM3'};
genes={'ATM', 'BCOR', 'CHD2', 'DDX3X', 'EGR2', 'FBXW7', 'HIST1H1E', 'ITPKB', 'KRAS', 'MED12', 'MYD88', 'NOTCH1', 'NRAS', 'POT1', 'RIPK1', 'SAMHD1', 'SF3B1', 'TP53', 'XPO1', 'ZMYM3', 'KLF4', 'RAB22A', 'SLC9A1', 'TOM1L2', 'WASF1', 'ACYP1', 'ATP11B', 'CLIP1', 'FANCE', 'HLAB', 'ICMT', 'JKAMP', 'MS4A4A', 'MUC5AC', 'NASP', 'NOL11', 'SAMHD1', 'SAMSN1', 'SGPP1', 'SPON1', 'SWSAP1', 'ZSWIM1', 'C3orf43', 'GAPDH', 'GNAS'};

%M=reorder_struct(M,~ismember(M.Patient,'J341'));

samples.names=unique(M.Patient);
%samples=reorder_struct(samples,~ismember(samples.names,'J341'));

N=length(samples.names);
NV=length(genes);



for v=1:NV
    for s=1:N
         k=ismember(M.Patient,samples.names{s});
         if ismember(genes{v},M.Gene(k));
             samples.(genes{v})(s,1)=1;
         else
             samples.(genes{v})(s,1)=0;
         end
    end
end



fields=fieldnames(samples);
samples=sort_struct(samples,{fields{2:end}},(ones(length({fields{2:end}}))*-1));
%define colors
relapse=[215/255 48/255 39/255];
diag=[39/255 48/255 215/255];

%add layers until you specify hold off
hold on


%N=100%length(unique(G.sample));
%NV=20%slength(si);
axis([-1.5 N+2 -2 NV+8]);
set(gca,'visible','off');
%setting dims of boxs 
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;

%double is to make string numeric
M.JustInRelapse=str2double(M.JustInRelapse);
for v=1:NV
    for s=1:N
     %these correspond to the shapes of the panels in this case rectangles   
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-v;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.95 .95 .95],'edgecolor',[.95 .95 .95])
        
    end
end

for v=1:NV
    for s=1:N
        % loop through the variants and samples placing patches for
        % mutations in either the germline half or somatic half based on
        % missense and nonsense
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];%
     j=NV-v;
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
     k=ismember(M.Patient,samples.names{s});
    if ismember(genes{v},M.Gene(k));
       
        l=find(ismember(M.Gene,genes{v})&ismember(M.Patient,samples.names{s}));
        if M.JustInRelapse(l)==1
        patch(x1,y1,[0 0 0],'edgecolor',[0 0 0],'EdgeAlpha',0)
        else
            patch(x1,y1,relapse,'edgecolor',relapse,'EdgeAlpha',0)
        end
    end
    end
end

for v=1:length(genes)
    j=NV-v;
    text (.5,j-.5,genes{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end