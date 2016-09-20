%BcellSignaling RNA expression plots
Bcell=load_struct('/Users/amaro/Downloads/BcellExpression_datamatrix.txt');
fields=fieldnames(Bcell);
for i=2:length(fields)
Bcell.(fields{i})=str2double(Bcell.(fields{i}));
end

%% patient 22
figure()
hold on
plot(Bcell.cll_xbam_tumor1,Bcell.cll_xbam_tumor2,'r.')
corrcoef(Bcell.cll_xbam_tumor1,Bcell.cll_xbam_tumor2)
text(100,500,'R=.96')
title('RNA exp of Bcell Signaling in patient 22')
xlabel('RPKM Tumor1')
ylabel('RPKM Tumor2')
h=refline(1,0);
set(h,'color','k','LineStyle','- -');
patient22_foldchange=Bcell.cll_xbam_tumor2./Bcell.cll_xbam_tumor1;
patient22_foldchange(isinf(patient22_foldchange))=NaN;
Bcell.Hugo(patient22_foldchange<.5)
%patient22 fold change: timeptB/timeptA
 %   'IKBKB'  %fold change:  0.490000000000000 chr8
  %  'PIK3R5' %fold change: 0.383259911894273 chr17

%% patient12
figure()
hold on
plot(Bcell.CLLMDAC0012TumorSM3VIED,Bcell.CLLMDAC0012TumorSM3VIEE,'r.')
corrcoef(Bcell.CLLMDAC0012TumorSM3VIED,Bcell.CLLMDAC0012TumorSM3VIEE);  
text(100,2000,'R=.99')
title('RNA exp of Bcell Signaling in patient 12')
xlabel('RPKM VIED')
ylabel('RPKM VIEE')
h=refline(1,0);
set(h,'color','k','LineStyle','- -');
patient12_foldchange=Bcell.CLLMDAC0012TumorSM3VIEE./Bcell.CLLMDAC0012TumorSM3VIED;
patient12_foldchange(isinf(patient12_foldchange))=NaN;
Bcell.Hugo(patient12_foldchange<.5)
%patient12  fold change: timeptB/timeptA
 %   'JUN' %fold change: 0.254883720930233 chr1
  %  'NFKBIA' %fold change: .45  chr14 
   
%% patient11
figure()
hold on
plot(Bcell.CLLMDAC0011TumorSM3VIEB,Bcell.CLLMDAC0011TumorSM4M92J,'r.')
a=corrcoef(Bcell.CLLMDAC0011TumorSM3VIEB,Bcell.CLLMDAC0011TumorSM4M92J);
text(100,350,num2str(a(2)));
title('RNA exp of Bcell Signaling in patient 11')
xlabel('RPKM VIEB')
ylabel('RPKM 4M92J')
h=refline(1,0);
set(h,'color','k','LineStyle','- -');
patient11_foldchange=Bcell.CLLMDAC0011TumorSM4M92J./Bcell.CLLMDAC0011TumorSM3VIEB;
%patient 11  fold change: timeptB/timeptA
%'PIK3R3' fold change: .15
 %   'JUN' fold change: .05
  %  'FOS' fold change: .002
   % 'LILRB3' fold change: .19
    %'RAC2' fold change.44


