files={'brain-SH0034-TP-NT-SM-5EV2E-SM-5EV2J_ABS_MAF.txt','brain-SH0420-TP-NT-SM-5EV2D-SM-5EV2I_ABS_MAF.txt','brain-SH0546-TP-NB-SM-4P7I2-SM-4VH4L_ABS_MAF.txt','brain-SH0622-TP-NT-SM-5EV2G-SM-5EV2L_ABS_MAF.txt',...
    'brain-SH0829-TP-NT-SM-5EV2F-SM-5EV2K_ABS_MAF.txt'...
,'brain-SH1537-TP-NB-SM-4P7HX-SM-4VH4H_ABS_MAF.txt','brain-SH3000-TP-NB-SM-4P7HY-SM-4VH4I_ABS_MAF.txt','brain-SH3979-TP-NB-SM-4P7HZ-SM-4VH4J_ABS_MAF.txt','brain-SH4885-TP-NT-SM-5EV2H-SM-5EV2M_ABS_MAF.txt','brain-SH85592-TP-NB-SM-4P7I1-SM-4VH4K_ABS_MAF.txt'}


files=files';
i=1;
M=load_struct(strcat('/Users/amaro/Downloads/HemangiomaABSOLUTE_Results/',files{i}));

for i=2:length(files)
    m=load_struct(strcat('/Users/amaro/Downloads/HemangiomaABSOLUTE_Results/',files{i}));
    M=mergeStruct(m,M);

end
M=rmfield(M,'N');
save_struct(M,strcat('/Users/amaro/Downloads/HemangiomaABSOLUTE_Results/CombinedMaf.tsv'));


for i=1:length(samples)
k=ismember(M.sample,samples{i});
pur(i,1)=unique(M.purity(k));
end

