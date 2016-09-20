% script to generate minimal sample info file from D

dbstop if error;
set_verbose_level(30);

% working directories
tcgascape_dir = '/xchip/gistic/tcgascape/';    % all runs
run_dir = [tcgascape_dir 'tcgascape_120416/']; % this run
igv_dir = [run_dir 'igvfiles/'];               % flat files

D = load_D([run_dir 'Dmaster.mat']);
cancer_namemap = tcga_cancer_types();
xtypes = values(cancer_namemap,{D.sis.gcmtype});
SI = struct('Array',D.sdesc','Disease',xtypes);
infofilewrite([igv_dir 'sample_info.txt'],SI);
