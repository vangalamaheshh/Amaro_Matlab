M=load_table('~/Documents/JMML_maf_for_suppl_fig3.txt');
M=load_table('~/Documents/JMML-UPN2531_1-TP-NB-SM-5861O-SM-5861N_Vs_JMML-UPN2531_1-TR-NB-SM-5861P-SM-5861N__pairwise_CCF_DP.maf.annotated.txt');
M.CCF_tot=M.CCF1+M.CCF2;
drivers={'SH2B3','PTPN11'};
M.isdriver=ismember(M.locus,drivers);
M=sort_struct(M,{'CCF_tot','isdriver'},[-1,-1]);

NV=slength(M);
N=2;
figure()
hold on
load('grncolormap.mat')

dx_cn=0.4; dy_cn=0.4;
for v=1:NV
    for s=1:N
        
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
   
        
        patch(x1,y1,grncolor(round((255*M.(strcat('CCF',num2str(s)))(v))+1),:)./255,'edgecolor',[0 0 0]);

   
    
        
    end
end

for v=1:NV
    s=-1;
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    text((s-dx_cn+s+dx_cn)/2,(j-dy_cn+j+dy_cn)/2,strcat(M.locus{v},'_',M.Protein_Change{v}),'HorizontalAlignment','center','FontSize',20);
    
end

set(gca,'visible','off');


%127,201,127