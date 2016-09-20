function C2=run_glad(prefix,C,rewrite_data,params)
%% FIXME: Why the heck can't we run this code unless we copy the
%getGLADparam.R into our own home directory?  When we change the
%~/R/getGLADparam.R to point to Gaddy's home directory, it can't read the
%file.  Unix read/write permissions do not appear to explain the
%problem.  Temporary but unacceptable fix is to copy R/getGLADparam.R into
%current user's home directory.
%
%           Revisions:
%               12 Dec 07 -- Fixed above problem (jdobson)
% ---
% $Id$
% $Date: 2007-09-18 13:20:16 -0400 (Tue, 18 Sep 2007) $
% $LastChangedBy: rameen $
% $Rev$

j
fname=[prefix '.dat'];
if exist(fname,'file')
  if ~exist('rewrite_data','var') || rewrite_data
    write_as_dchip(fname,C);
  end
else
    write_as_dchip(fname,C);
end

f=fopen([prefix '.R'],'w');
if exist('params','var') && ~isempty(params)
  fprintf(f,['source("~/CancerGenomeAnalysis/trunk/R/getGLADparam.R")' newline]);
  st=[];
  if isfield(params,'qlambda')
    st=[st ',ql=' num2str(params.qlambda)];
  end
  if isfield(params,'bandwidth')
    st=[st ',bw=' num2str(params.bandwidth)];
  end
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
  fprintf(f,['getGLAD("' prefix '.dat","' prefix '.seg.dat"' st ')' newline]);
else
  fprintf(f,['source("~/CancerGenomeAnalysis/trunk/R/getGLAD.R")' newline]);
  fprintf(f,['getGLAD("' prefix '.dat","' prefix '.seg.dat")' newline]);
end
fprintf(f,['q()' newline]);
fclose(f);
