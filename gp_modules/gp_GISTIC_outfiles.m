function gp_GISTIC_outfiles(varargin)
% gp_GISTIC_outfiles  -ui path to input.mat -GISTIC
% GISTICout.mat -p param file -[-basedir base_directory] 
% -lout all_lesions.txt  -sout scores.txt -add_vals 1expandsall_lesion
%  -[ -name my_struct ] -[-refgene_file  hg17_20070131.mat] -[-annot_file
%  geneannotations] -[-mark annotation_marker]
% ---
% $Id$
% $Date: 2007-09-18 13:20:16 -0400 (Tue, 18 Sep 2007) $
% $LastChangedBy: rameen $
% $Rev$
  
addpath ~/CancerGenomeAnalysis/trunk/matlab/
addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp

a=handle_args({'ui','GISTIC','p','basedir', 'lout','sout','add_vals', 'name', 'refgene_file', 'annot_file', 'mark'},varargin);


GISTIC_outfiles(a);


% for compilation:
% cd ~/matlab/gp_modules
% mcc -m gp_GISTIC_outfiles
