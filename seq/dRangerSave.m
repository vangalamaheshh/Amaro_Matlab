function dRangerSave(sample,T,N)

x = T;
file = ['/xchip/tcga_scratch/lawrence/' sample '/tumor_dR/all.weird.joined.mat'];
save(file,'x','-v7.3');

x = N;
file = ['/xchip/tcga_scratch/lawrence/' sample '/normal_dR/all.weird.joined.mat'];
save(file,'x','-v7.3');

