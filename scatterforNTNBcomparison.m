NtNB=load_table('/Users/amaro/Documents/deTiN_Gaddy_LabMeeting/TableForNTNBcomparisonalmostfilledin.txt')
% 
%     'BLCA' 255?170?170
%     'BRCA' 255?209?170
%     'COAD' 102?153?153
%     'HNSC' 136?204?136
%     'KIRC' 170?57?57
%     'KIRP' 170?108?57
%     'LIHC'
%     'LUAD'
%     'LUSC'
%     'PRAD'
%     'STAD'
%     'THCA'
% 
%    shade 0 = #F62337 = rgb(246, 35, 55) = rgba(246, 35, 55,1) = rgb0(0.965,0.137,0.216)
%    shade 1 = #F9D1D5 = rgb(249,209,213) = rgba(249,209,213,1) = rgb0(0.976,0.82,0.835)
%    shade 2 = #F76976 = rgb(247,105,118) = rgba(247,105,118,1) = rgb0(0.969,0.412,0.463)
%    shade 3 = #D60014 = rgb(214,  0, 20) = rgba(214,  0, 20,1) = rgb0(0.839,0,0.078)
%    shade 4 = #AA0010 = rgb(170,  0, 16) = rgba(170,  0, 16,1) = rgb0(0.667,0,0.063)
% 
% *** Secondary color (1):
% 
%    shade 0 = #FF7E24 = rgb(255,126, 36) = rgba(255,126, 36,1) = rgb0(1,0.494,0.141)
%    shade 1 = #FFE7D6 = rgb(255,231,214) = rgba(255,231,214,1) = rgb0(1,0.906,0.839)
%    shade 2 = #FFA96C = rgb(255,169,108) = rgba(255,169,108,1) = rgb0(1,0.663,0.424)
%    shade 3 = #DD5C00 = rgb(221, 92,  0) = rgba(221, 92,  0,1) = rgb0(0.867,0.361,0)
%    shade 4 = #AF4900 = rgb(175, 73,  0) = rgba(175, 73,  0,1) = rgb0(0.686,0.286,0)
% 
% *** Secondary color (2):
% 
%    shade 0 = #16A08D = rgb( 22,160,141) = rgba( 22,160,141,1) = rgb0(0.086,0.627,0.553)
%    shade 1 = #A0BEBA = rgb(160,190,186) = rgba(160,190,186,1) = rgb0(0.627,0.745,0.729)
%    shade 2 = #46A497 = rgb( 70,164,151) = rgba( 70,164,151,1) = rgb0(0.275,0.643,0.592)
%    shade 3 = #008B78 = rgb(  0,139,120) = rgba(  0,139,120,1) = rgb0(0,0.545,0.471)
%    shade 4 = #006E5F = rgb(  0,110, 95) = rgba(  0,110, 95,1) = rgb0(0,0.431,0.373)
% 
% *** Complement color:
% 
%    shade 0 = #3ED41E = rgb( 62,212, 30) = rgba( 62,212, 30,1) = rgb0(0.243,0.831,0.118)
%    shade 1 = #C4E2BE = rgb(196,226,190) = rgba(196,226,190,1) = rgb0(0.769,0.886,0.745)
%    shade 2 = #70D75B = rgb(112,215, 91) = rgba(112,215, 91,1) = rgb0(0.439,0.843,0.357)
%    shade 3 = #20B900 = rgb( 32,185,  0) = rgba( 32,185,  0,1) = rgb0(0.125,0.725,0)
%    shade 4 = #199200 = rgb( 25,146,  0) = rgba( 25,146,  0,1) = rgb0(0.098,0.573,0)   
   
tumor_types=unique(NtNB.Tumor_Type);
%color_tumor=[0.667,0.224,0.224;0.333,0,0;0.667,0.424,0.224;0.333,0.153,0;0.133,0.4,0.4;0,0.2,0.2;0.176,0.533,0.176;0,0.267,0;0.251,0.498,0.498;0.333,0.667,0.333;0.831,0.416,0.416;0.831,0.604,0.416];
color_tumor=[0.965,0.137,0.216;0.969,0.412,0.463;0.667,0,0.063;1,0.494,0.141;1,0.663,0.424;0.686,0.286,0;0.086,0.627,0.553;0.275,0.643,0.592;0,0.431,0.373;0.243,0.831,0.118;0.439,0.843,0.357;0.098,0.573,0];

NtNB=sort_struct(NtNB,'tumor_in_normal_estimate',-1);
kNT=ismember(NtNB.comparison,'TP-NT');
kNB=ismember(NtNB.comparison,'TP-NB');

color_array=zeros(slength(NtNB),3);


%kNT
for i=1:length(tumor_types)
color_array(ismember(NtNB.Tumor_Type,tumor_types{i}),1)=color_tumor(i,1);
color_array(ismember(NtNB.Tumor_Type,tumor_types{i}),2)=color_tumor(i,2);
color_array(ismember(NtNB.Tumor_Type,tumor_types{i}),3)=color_tumor(i,3);
end

stepNT=.2/(sum(kNT)-1);
stepNB=.2/(sum(kNB)-1);

scatter((1-.1):stepNB:((1+.1)),NtNB.tumor_in_normal_estimate(kNB),25,color_array(kNB,:),'filled')
hold on
scatter((2-.1):stepNT:((2+.1)),NtNB.tumor_in_normal_estimate(kNT),25,color_array(kNT,:),'filled')
set(gca,'XTick',[1:2],'XTickLabel',{'Normal Blood','Normal Tissue'})
ylim([0 .5])
xlim([.5,2.5])

for i=1:length(tumor_types)
text (.75,.4-(.02*i),tumor_types{i},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)

end
dy_cn=.005;
dx_cn=.01;
j=.39;
s=.77;
plot([.25;2.75],[.02;.02],'k--')
for i=1:length(tumor_types)
    x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    patch(x1,y1,color_tumor(i,:),'edgecolor',color_tumor(i,:))
    j=j-.02;
end
