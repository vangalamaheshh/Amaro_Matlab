%TSPSSignaling RNA expression plots
TSPS=load_struct('/Users/amaro/Downloads/TSPSmatrix.txt');
fields=fieldnames(TSPS);
for i=2:length(fields)
TSPS.(fields{i})=str2double(TSPS.(fields{i}));
end

%% patient 22
figure()
hold on
plot(TSPS.cll_xbam_tumor1,TSPS.cll_xbam_tumor2,'r.')
corrcoef(TSPS.cll_xbam_tumor1,TSPS.cll_xbam_tumor2)
text(100,500,'R=.96')
title('RNA exp of TSPS Signaling in patient 22')
xlabel('RPKM Tumor1')
ylabel('RPKM Tumor2')
h=refline(1,0);
set(h,'color','k','LineStyle','- -');
patient22_foldchange=TSPS.cll_xbam_tumor2./TSPS.cll_xbam_tumor1;
patient22_foldchange(isinf(patient22_foldchange))=NaN;
TSPS.Hugo(patient22_foldchange<.5)

%% patient12
figure()
hold on
plot(TSPS.CLLMDAC0012TumorSM3VIED,TSPS.CLLMDAC0012TumorSM3VIEE,'r.')
corrcoef(TSPS.CLLMDAC0012TumorSM3VIED,TSPS.CLLMDAC0012TumorSM3VIEE);  
text(100,2000,'R=.99')
title('RNA exp of TSPS Signaling in patient 12')
xlabel('RPKM VIED')
ylabel('RPKM VIEE')
h=refline(1,0);
set(h,'color','k','LineStyle','- -');
patient12_foldchange=TSPS.CLLMDAC0012TumorSM3VIEE./TSPS.CLLMDAC0012TumorSM3VIED;
patient12_foldchange(isinf(patient12_foldchange))=NaN;
TSPS.Hugo(patient12_foldchange<.5)
   
%% patient11
figure()
hold on
plot(TSPS.CLLMDAC0011TumorSM3VIEB,TSPS.CLLMDAC0011TumorSM4M92J,'r.')
a=corrcoef(TSPS.CLLMDAC0011TumorSM3VIEB,TSPS.CLLMDAC0011TumorSM4M92J);
text(100,350,num2str(a(2)));
title('RNA exp of TSPS Signaling in patient 11')
xlabel('RPKM VIEB')
ylabel('RPKM 4M92J')
h=refline(1,0);
set(h,'color','k','LineStyle','- -');
patient11_foldchange=TSPS.CLLMDAC0011TumorSM4M92J./TSPS.CLLMDAC0011TumorSM3VIEB;
patient11_foldchange(isinf(patient12_foldchange))=NaN;
TSPS.Hugo(patient11_foldchange<.5)


