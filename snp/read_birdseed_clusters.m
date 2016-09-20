function B=read_birdseed_clusters(fname)

form='%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';

f=fopen(fname);
d=textscan(f,form,'bufSize',1e7,'delimiter',' ;');
fclose(f);

B.marker=d{1};
B.dat=cat(2,d{2:end});
B.sdesc={'AA_mean_a','AA_mean_b','AA_covar_aa','AA_covar_ab','AA_covar_bb','AA_num_observations',...
         'AB_mean_a','AB_mean_b','AB_covar_aa','AB_covar_ab','AB_covar_bb','AB_num_observations',...
         'BB_mean_a','BB_mean_b','BB_covar_aa','BB_covar_ab','BB_covar_bb','BB_num_observations'};
clear d
