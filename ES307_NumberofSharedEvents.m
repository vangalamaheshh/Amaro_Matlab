EAC1=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TP-NT-SM-4P8W1-SM-4QZMM.pon_filtered.txt');
ind=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TP-NT-SM-4P8W1-SM-4QZMM.indel.capture.maf.annotated');
EAC1=mergeStruct(EAC1,ind);
EAC1.key=strcat(EAC1.Chromosome,EAC1.Start_position);
EAC3=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TP-NT-SM-4QZMO-SM-4QZMM.pon_filtered.txt');
ind=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TP-NT-SM-4QZMO-SM-4QZMM.indel.capture.maf.annotated');
EAC3=mergeStruct(EAC3,ind);
HGD1=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TP-NT-SM-4P8VZ-SM-4QZMM.pon_filtered.txt');
ind=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TP-NT-SM-4P8VZ-SM-4QZMM.indel.capture.maf.annotated');
HGD1=mergeStruct(HGD1,ind);
MET1=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TM-NT-SM-4P8W3-SM-4QZMM.pon_filtered.txt');
ind=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ES_3-07_Gaddy/ES-3_07-TM-NT-SM-4P8W3-SM-4QZMM.indel.capture.maf.annotated');
MET1=mergeStruct(MET1,ind);
EAC3.key=strcat(EAC3.Chromosome,EAC3.Start_position);
MET1.key=strcat(MET1.Chromosome,MET1.Start_position);
HGD1.key=strcat(HGD1.Chromosome,HGD1.Start_position);
