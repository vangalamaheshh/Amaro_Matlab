function maf_split(mdir)
x=dir(mdir);
for i=3:length(x)
m=load_struct(strcat(mdir,x(i).name));
sm=reorder_struct(m,ismember(m.Variant_Type,'SNP')|ismember(m.Variant_Type,'DNP'));
im=reorder_struct(m,ismember(m.Variant_Type,'DEL')|ismember(m.Variant_Type,'INS'));
s=regexp(x(i).name,'.maf','split');
base=s{1};

save_struct(sm,strcat(mdir,base,'.Split_SNVs.maf'));
save_struct(im,strcat(mdir,base,'.Split_Indels.maf'));
end


end

