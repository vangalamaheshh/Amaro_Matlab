function [M,C,P,batches,snps,mn]=gp_batch_correct(infile,sample_info,outfile,min_sz,bonf_pv_thresh,absolute_pv)

if ~exist('min_sz','var')
  min_sz=5;
end

if ~exist('bonf_pv_thresh','var')
  bonf_pv_thresh=0.05;
end

if ~exist('absolute_pv','var')
  absolute_pv=0.001;
end
  
M=read_modelled_data_file(infile,-1,1,1);
SI=read_sample_info_file(sample_info);
M=add_samples_info(M,SI,'name');
% take log2 if needed
if median(M.dat(:))>0.5
  disp('Taking Log2 of the data after thresholding at 1');
  M.dat(M.dat<1)=1;
  M.dat=log2(M.dat)-1; % not really needed
end
disp('Correcting batch effect');
[C,P,batches,snps,mn]=correct_batch_effect(M,min_sz,bonf_pv_thresh,absolute_pv);
disp('Taking 2^ of the data');
C.dat=2.^(C.dat+1); % back to copy number

write_as_dchip(outfile,CL,1);


