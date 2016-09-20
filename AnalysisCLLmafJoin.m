for i=1:length(I)
k=find(ismember(ABS.individual,I{i}));
for j=1:length(k)
if j==1
m=load_struct(ABS.ABS_Maf{k(j)});
M.Gene=m.Hugo_Symbol;
M.Start_position=m.Start_position;
patient=m.Tumor_Sample_UUID{1};
patient=tr(patient,'-','_');
M.(strcat(patient))=m.ccf_hat;
else
    m=load_struct(ABS.ABS_Maf{k(j)});
    patient=m.Tumor_Sample_UUID{1};
    patient=tr(patient,'-','_');
    M.(strcat(patient))=m.ccf_hat;
end
end
save_struct(M,strcat(I{i},'_ABSMAF.txt'));
end