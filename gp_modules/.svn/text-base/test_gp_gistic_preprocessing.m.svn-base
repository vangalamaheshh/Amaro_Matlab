if ispc
  addpath('K:\CancerGenomeAnalysis\trunk\matlab');
  rmpath_regexp('CancerGenomeAnalysis');
  addpath('K:\CancerGenomeAnalysis\trunk\matlab');
  addpath('K:\CancerGenomeAnalysis\trunk\matlab\snp');
  addpath('K:\CancerGenomeAnalysis\trunk\matlab\gp_modules');
  addpath('K:\CancerGenomeAnalysis\trunk\matlab\HDF5');
  base_dir = 'N:\GPTests\GISTICPreprocessing\';  %b
else
  addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab');
  rmpath_regexp('CancerGenomeAnalysis');
  addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab');
  addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab/snp');
  addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab/HDF5');
  addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab/gp_modules');
  base_dir = '/home/radon00/cangen/GPTests/GISTICPreprocessing/';  %b
end
array_list_file = [base_dir 'array_list.txt'];  %al
sample_info_file = [base_dir 'glioma_si_070424.xls' ]; %si
bc = '1';  %bc
snp_infiles = ...
    [ base_dir '100H_Glioma_070921.snp ' base_dir '100X_Glioma_070921.snp' ];
    %i
snp_outfile =  [ base_dir 'gliomasnps2.snp']; %o
save_mat = '1';  %sm
histqctumors = '1';
show_hist = '0';
save_raw = '0';

set_verbose_level(30);

if ispc
  disp(['cd N:\GPTests\GISTICPreprocessing\']);
  disp(['K:\CancerGenomeAnalysis\trunk\matlab\gp_modules\gp_gistic_preprocessing.exe ' ...
        '-al ' 'array_list.txt ' '-si ' 'glioma_si_070424.xls ' ...
        '-b ' 'N:\GPTests\GISTICPreprocessing\ ' '-bc ' '1 ' ...
        '-i ' '100H_Glioma_070921.snp 100X_Glioma_070921.snp ' '-o ' 'gliomasnps2.snp ' '-sm ' '1 '  ...
        '-ht ' '1 ' '-sh ' '0 ' '-svr ' '0 ']);

else
  gp_gistic_preprocessing('-al',array_list_file,'-si',sample_info_file,...
    '-b',base_dir,'-bc',bc,'-i',snp_infiles,'-o',snp_outfile,'-sm',save_mat,...
    '-ht', histqctumors,'-sh',show_hist,'-svr',save_raw)
end

%To do:  add check of file/directory read/write permissions early in
%gp_gistic preprocessing.

