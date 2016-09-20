function gp_GISTIC_armfigures(varargin)
% gp_GISTIC_outfiles  -ui path to input.mat -GISTIC
% GISTICout.mat -p param file -[-basedir base_directory] -n number of samples
% ---
% $Id$
% $Date: 2007-09-17 15:18:13 -0400 (Mon, 17 Sep 2007) $
% $LastChangedBy: barbara$
% $Rev$

addpath ~/CancerGenomeAnalysis/trunk/matlab/
addpath ~/CancerGenomeAnalysis/trunk/matlab/gp_modules
addpath ~/CancerGenomeAnalysis/trunk/matlab/snp


a=handle_args({'ui','GISTIC','p','basedir'},varargin);


GISTIC_armfigures(a);


% for compilation:
% cd ~/matlab/gp_modules
% mcc -m gp_GISTIC_outfiles
