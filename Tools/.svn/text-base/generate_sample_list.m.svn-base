cd /xchip/gistic/GCM/preprocessing_doubletan_nobc_pinv_080328/output

D = load_DWP('D.mat');  %load datastruct in "write-protect" mode (avoids having to transfer hdf5 files)

names = D.sdesc;
normalsbool = D.supdat(strmatch('N',D.supacc,'exact'),:);
sample_number = 1:length(names);

fid = fopen('SampleOrders.txt','w');

for k = 1:length(names)
fprintf(fid,[names{k} '\t\t%d\t%d\n'],normalsbool(k),sample_number(k));
end

fclose(fid);