x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0012-CCFdist.txt',' ');

for i=1:size(x,2)
    x{i}=str2double(x{i});
    
    x1=reshape(x{i},101,101);
    TP1.(strcat('MDAC_0012_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0012_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end

x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0022_1_2-CCFdist.txt',' ');




for i=1:size(x,2)
    x1=str2double(x{i});
    x1=reshape(x1,101,101);
    TP1.(strcat('MDAC_0022_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0022_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end


x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0022_3_4-CCFdist.txt',' ');
for i=1:size(x,2)
    x1=str2double(x{i});
    x1=reshape(x1,101,101);
    TP1.(strcat('MDAC_0022_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0022_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end


x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0011_K_J-CCFdist.txt',' ');
for i=1:size(x,2)
    x1=str2double(x{i});
    x1=reshape(x1,101,101);
    TP1.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end


x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0011_K_C-CCFdist.txt',' ');
for i=1:size(x,2)
    x1=str2double(x{i});
    x1=reshape(x1,101,101);
    TP1.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end




x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0011_J_K-CCFdist.txt',' ');
for i=1:size(x,2)
    x1=str2double(x{i});
    x1=reshape(x1,101,101);
    TP1.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end




x=read_dlm_file('/Users/amaro/Downloads/Rdata_forCLLRush/MDAC-0011_K_C-CCFdist.txt',' ');
for i=1:size(x,2)
    x1=str2double(x{i});
    x1=reshape(x1,101,101);
    TP1.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,2);
    [v,l]=max(sum(x1,2));
    clone_locs(i,1)=l;
    TP2.(strcat('MDAC_0011_clone',num2str(i)))=sum(x1,1)';
    [v,l]=max(sum(x1,1));
     clone_locs(i,2)=l;
end

