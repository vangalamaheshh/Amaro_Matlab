m=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRLX-SM-5MRLZ_Vs_ES-7_10-TP-NT-SM-5MRLU-SM-5MRLZ__pairwise_CCF_DP.maf.annotated.txt');
CCFDP_m=m;
m2=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRLU-SM-5MRLZ.pon_filtered.txt');
m1=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRLX-SM-5MRLZ.pon_filtered.txt');
m3=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRLW-SM-5MRLZ.pon_filtered.txt');
m4=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRLV-SM-5MRLZ.pon_filtered.txt');
m5=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRM1-SM-5MRLZ.pon_filtered.txt');
m6=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRM2-SM-5MRLZ.pon_filtered.txt');
m7=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRM3-SM-5MRLZ.pon_filtered.txt');
m8=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRM4-SM-5MRLZ.pon_filtered.txt');
m9=load_struct('/Users/amaro/Downloads/ES-7_10-TP-NT-SM-5MRLY-SM-5MRLZ.pon_filtered.txt');
CCFDP_m.key=strcat(CCFDP_m.Chr,CCFDP_m.pos);
m1.key=strcat(m1.Chromosome,m1.Start_position);
m2.key=strcat(m2.Chromosome,m2.Start_position);

M=mergeStruct(m1,m2);
M=mergeStruct(M,m3);
M=mergeStruct(M,m4);
M=mergeStruct(M,m5);

M=mergeStruct(M,m6);
M=mergeStruct(M,m7);

M=mergeStruct(M,m8);
M=mergeStruct(M,m9);
M.key=strcat(M.Chromosome,M.Start_position);

violate_keys=CCFDP_m.key(ismember(CCFDP_m.var_cluster_classes,'S_S'));

