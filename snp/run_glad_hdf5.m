function run_glad_hdf5(chrposfile,prefix,hdf5file,rowidx,params)
% RUN_GLAD_HDF5 write R script to invoke getGLAD_hdf5 which load sample
% data directly from hdf5.
%
%   C2 = RUN_GLAD_HDF5(PREFIX,HDF5FILE,ROWIDX,PARAMS);
%
%           File added 15 Feb 08  - jdobson

% ---
% $Id$
% $Date: 2008-04-25 11:16:44 -0400 (Fri, 25 Apr 2008) $
% $LastChangedBy: jdobson $
% $Rev$



f=fopen([prefix '.R'],'w');
if exist('params','var') && ~isempty(params)
  fprintf(f,['source("~/CancerGenomeAnalysis/trunk/R/getGLADparam_hdf5.R")' newline]);
  st=[];

  if isfield(params,'lambdabreak')
    st=[st ',lb=' num2str(params.lambdabreak)];
  end
  if isfield(params,'lambdacluster')
    st=[st ',lc=' num2str(params.lambdacluster)];
  end
  if isfield(params,'lambdaclusterGen')
    st=[st ',lcg=' num2str(params.lambdaclusterGen)];
  end
  if isfield(params,'param_d')
    st=[st ',pr=c(d=' num2str(params.param_d) ')'];
  end 
  if isfield(params,'qlambda')
      st=[st ',ql=' num2str(params.qlambda)];
  end
  if isfield(params,'bandwidth')
      st=[st ',bw=' num2str(params.bandwidth)];
  end
  fprintf(f,['getGLADparam_hdf5("' chrposfile '","' prefix '.seg.dat", "' ...
      hdf5file '", "dat", ' num2str(rowidx)  st ')' newline]);
else
  fprintf(f,['source("~/CancerGenomeAnalysis/trunk/R/getGLAD_hdf5.R")' newline]);
  fprintf(f,['getGLAD_hdf5("' chrposfile '","' prefix '.seg.dat", "' ...
      hdf5file '" ,"dat",' num2str(rowidx) ')' newline]);  %Note: st is not passed in
end
fprintf(f,['q()' newline]);
fclose(f);
