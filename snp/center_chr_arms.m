function [dat, MEDS] = center_chr_arms(C, cyto)

disp('adding cyto info...');
C=add_cyto(C, cyto);
C.chrarmn=C.armn+C.chrn*2-2;
dat=C.dat;
s=size(C.dat);
meds=zeros(48,s(2));
MEDS=zeros(1, s(2));
disp('finding arm medians')
for i = 1:max(C.chrarmn)
    disp(i);
    meds((i),:)=nanmedian(dat(find(C.chrarmn==i),:));
    s= size(dat(find(C.chrarmn==i),:));
    medsubtractor=repmat(meds((i),:), s(1),1);
    dat(find(C.chrarmn==i),:)=dat(find(C.chrarmn==i),:) - medsubtractor;
    MEDS(find(C.chrarmn==i),:)=medsubtractor;
end

