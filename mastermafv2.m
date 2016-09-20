files=dir('/Users/amaro/Downloads/MasterMAFv2/')
files.name
files(1)=[];
files(end)=[];
files.name
m=load_struct(strcat('/Users/amaro/Downloads/MasterMAFv2/',files(1).name));
MM=m
for i=2:length(files)
m=load_struct(strcat('/Users/amaro/Downloads/MasterMAFv2/',files(i).name));
MM=mergeStruct(MM,m);
end
k=cellfun(@isempty,MM.Tumor_Sample_UUID);
sum(k)
MM.Tumor_Sample_UUID(k)=MM.sample(k);
k=cellfun(@isempty,MM.sample);
MM.sample(k)=MM.Tumor_Sample_UUID(k);
MM=rmfield(MM,'N');
MM=reorder_struct(MM,~ismember(MM.sample,'NA'));

MM.x1{836}=MM.x099{836};
MM.x1
MM.key=strcat(MM.sample,'-',MM.Start_position,'-',MM.Tumor_Seq_Allele2);
save_struct(MM,'/Users/amaro/Downloads/MasterMAFv2/MasterMafFile.txt')


files=dir('/Users/amaro/Downloads/MasterMAFv2/mafs_for_annotation');
files(1)=[];
files(1)=[];
m=load_struct(strcat('/Users/amaro/Downloads/MasterMAFv2/mafs_for_annotation/',files(1).name));
f_val=fieldnames(m);
f_val
f_val{end-22:end-1}
    m.key=strcat(m.Tumor_Sample_UUID,'-',m.Start_position,'-',m.Tumor_Seq_Allele2);

validation_fields={f_val{end-20:end-1}};

for i=1:11%validation mafs 
    
    m=load_struct(strcat('/Users/amaro/Downloads/MasterMAFv2/mafs_for_annotation/',files(i).name));
    m.key=strcat(m.Tumor_Sample_UUID,'-',m.Start_position,'-',m.Tumor_Seq_Allele2);
    
    for j=1:slength(m)
        k=find(ismember(MM.key,m.key{j}));
        
        if ~isempty(k)
            k=k(1);
            for z=1:length(validation_fields)
                if isfield(m,validation_fields{z})
                MM.(validation_fields{z}){k,1}=m.(validation_fields{z}){j,1};
                else
                     MM.(validation_fields{z}){k,1}='NA';
                end
            end
        end
    
    
    end
end
m=load_struct(strcat('/Users/amaro/Downloads/MasterMAFv2/mafs_for_annotation/',files(12).name));
f_val=fieldnames(m);
f_val
f_val{end-166:end-1}
    m.key=strcat(m.Tumor_Sample_UUID,'-',m.Start_position,'-',m.Tumor_Seq_Allele2);

validation_fields={f_val{end-166:end-1}};

k=ismember(MM.Matched_Norm_Sample_Barcode,'HSC-UTMDAC-0001-Normal-SM-4ET5L');
for z=1:length(validation_fields)
    MM.(validation_fields{z})(k)={'NA'};
end

for i=12:length(files)
    m=load_struct(strcat('/Users/amaro/Downloads/MasterMAFv2/mafs_for_annotation/',files(i).name));
        m.key=strcat(m.sample,'-',m.Start_position,'-',m.Tumor_Seq_Allele2);
        for j=1:slength(m)
        k=find(ismember(MM.key,m.key{j}));
        
        if ~isempty(k)
            k=k(1);
            for z=1:length(validation_fields)
                if isfield(m,validation_fields{z})
                MM.(validation_fields{z}){k,1}=m.(validation_fields{z}){j,1};
                else
                     MM.(validation_fields{z}){k,1}='NA';
                end
            end
        end
    
    
    end
end


