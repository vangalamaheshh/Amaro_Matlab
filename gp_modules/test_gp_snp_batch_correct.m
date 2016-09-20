addpath ~/matlab/snp
addpath ~/matlab/gp_modules


cd ~/projects/snp/pipeline


[M,C,P,hlarge,snps]=gp_snp_batch_correct('modelInput.txt','~/projects/snp/Rameen/data/sample info_051220.txt',...
                                         'test.out.txt'); %outfile,min_sz,bonf_pv_thresh,absolute_pv

gp_snp_batch_correct modelInput.txt ~/projects/snp/Rameen/data/sample\ info_051220.txt test.out.txt

