

 M=load_struct('/Volumes/xchip_cga_home/amaro/TumorInNormal/AML_2.0/AggregateAML.maf');
 sig_genes=load_struct('~/Documents/SigGeneListForAML.txt');
sig_genes=reorder_struct(sig_genes,~ismember(sig_genes.gene,{'ADAM30';'NHS';'TTK'})); %genes which are not expressed / did not validate 

% M=load_struct('/Volumes/xchip_cga_home/amaro/CLL/DeTiN_Figure/DeTiNMAF_PanCan.maf');%panCan Maf
% sig_genes=load_struct('~/Documents/CGC_PanCanGenes.txt');

figure()
M=reorder_struct(M,ismember(M.Hugo_Symbol,sig_genes.gene));
sig_genes=reorder_struct(sig_genes,ismember(sig_genes.gene,M.Hugo_Symbol));
samples.names=unique(M.Tumor_Sample_Barcode);

N=length(samples.names);
NV=length(sig_genes.gene);
for i=1:slength(sig_genes)
    sig_genes.sum(i,1)=sum(ismember(M.Hugo_Symbol,sig_genes.gene{i}));
end   
sig_genes=sort_struct(sig_genes,'sum',-1);
for v=1:NV
    for s=1:N
         k=ismember(M.Tumor_Sample_Barcode,samples.names{s});
         if ismember(sig_genes.gene{v},M.Hugo_Symbol(k));
             samples.(sig_genes.gene{v})(s,1)=1;
         else
             samples.(sig_genes.gene{v})(s,1)=0;
         end
        
    end
    sig_genes.sum(v,1)=sum(samples.(sig_genes.gene{v}));
end



fields=fieldnames(samples);
samples=sort_struct(samples,{fields{2:end}},(ones(length({fields{2:end}}))*-1));
deTiN=[215/255 48/255 39/255];
normal=[0 0 0];


hold on


%N=100%length(unique(G.sample));
%NV=20%slength(si);
axis([-1.5 N+2 -2 NV+8]);
set(gca,'visible','off');
%setting dims of boxs 
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.38; dy_mut=.38;

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
     k=ismember(M.Tumor_Sample_Barcode,samples.names{s});
    if ismember(sig_genes.gene{v},M.Hugo_Symbol(k));
       
        l=find(ismember(M.Hugo_Symbol,sig_genes.gene{v})&ismember(M.Tumor_Sample_Barcode,samples.names{s}));
        if isequal(M.i_failure_reasons(l),{''})||isequal(M.i_failure_reasons(l),{',PoN'})
        patch(x1,y1,normal,'edgecolor',normal,'EdgeAlpha',0)
        else
            patch(x1,y1,deTiN,'edgecolor',deTiN,'EdgeAlpha',0)
        end
       
    end
    
    end
end

for v=1:slength(sig_genes)
    j=NV-v;
    text (.5,j-.5,sig_genes.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end