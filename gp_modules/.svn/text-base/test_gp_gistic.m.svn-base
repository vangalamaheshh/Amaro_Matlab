rmpath_regexp('trunk');
addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab');
addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab/snp');
addpath('/home/radon00/cangen/CancerGenomeAnalysis/trunk/matlab/gp_modules');
set_verbose_level(30);

basedir = '/home/radon00/cangen/GPTests/GISTIC/';
%segfile = [ basedir 'gliomasnps_cbs.no_header.seg'];
%segfile = [ basedir 'gliomasnps.glad'];
segfile = [ basedir 'output.no_header_appandcut071210_2.seg' ];%the segmentation file from March run

snpfile = [ basedir 'gliomasnps.snp'];
refgene = [ basedir 'hg16_20070112.mat'];
cnvfile = [ basedir 'Xba_Hind_CVN.061218.txt' ];
array_list_file = [ basedir 'core_array_list_071208.no_cell_lines.txt'];
segheaderlines = 0;


gp_gistic('-b',basedir,'-sgfl',segfile,'-snpfile',snpfile,'-refgene',...
    refgene,'-alf',array_list_file,'-cnv',cnvfile,'-ext','.oldseg');

