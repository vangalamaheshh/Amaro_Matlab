function gp_gistic_plots(varargin)
%Calls gistic_plots to creates GISTIC figures for amplification & deletion based on pipeline data
%that is stored in CL21, {regs, stats, and cyto}=GISTIC_core output matlab files
% ---
% $Id$
% $Date: 2007-12-03 15:31:14 -0500 (Mon, 03 Dec 2007) $
% $LastChangedBy: rameen $
% $Rev$

  
%addpath ~zhou/Matlab/Copied_Files
%addpath ~zhou/Matlab/Copied_Files/snp
%addpath ~zhou/Matlab/Plot_Wrap

addpath /xchip/gistic/Code/GISTIC1.0/matlab/
addpath /xchip/gistic/Code/GISTIC1.0/matlab/gp_modules
addpath /xchip/gistic/Code/GISTIC1.0/matlab/snp


a=handle_args({'uiCL21','GISTIC','p','basedir', 'uiCyto', 'ext','qt','ull','wp','ns','abb','as','ds','aLabels','dLabels','aTopLabels','dTopLabels','gn'},varargin);

%at least one of them has a base_dir (output_dir from previous module)
if ~isempty(a.basedir)
    base_dir = a.basedir;
    params=[];
else
  if ~isempty(a.p)
    params=read_params_file(a.p);
    pidx=grep('^output_dir$',{params.param},1);
    if isempty(pidx)
      error('No previous output_dir in params file');
    else
      base_dir=params(pidx(1)).value;
    end
  else
    error('Need to provide either a base directory or parameter file!');
  end
end

Gcore=[base_dir a.GISTIC];
G=load(Gcore);
   
if isfield(G,'cyto')
  cyto=G.cyto;
end

if isfield(G,'regs')
  regs=G.regs;
end

if isfield(G,'ts')
  ts=G.ts;
end

if isfield(G,'pvs')
  pvs=G.pvs;
end

if isfield(G,'p')
  p=G.p;
end

if isfield(G,'q')
  q=G.q;
end

if isfield(G,'ads')
  ads=G.ads;
end

if isfield(G,'d')
  d=G.d;
end

if isempty(a.uiCyto)
  error('Missing input cyto file name');
end
infileCyto=a.uiCyto;

disp(['Reading input file: ' infileCyto]);
load(infileCyto);
xCyto=load(infileCyto);

if isempty(xCyto)
   error('No variables in input cyto file');
end
  
if isempty(a.uiCL21)
  error('Missing input CL21 file name');
end
infileCL21=a.uiCL21;

disp(['Reading input file: ' infileCL21]);
xCL21=load(infileCL21);
if ~exist('xCL21')
  error('No variables in input CL21 file');
end

if isempty(a.ull)
    a.ull = 1;
else
    a.ull = str2num(a.ull);
end

if isempty(a.qt)
    a.qt = 0.25;
else
    a.qt = str2num(a.qt);
end

if isempty(a.wp)
    a.wp = 1;
else
    a.wp = str2num(a.wp);
end

if isempty(a.ns)
    a.ns = 1;
else
    a.ns = str2num(a.ns);
end

if isempty(a.abb)
    a.abb = 0;
else
    a.abb = str2num(a.abb);
end

if ~isempty(a.as) && ~isempty(a.ds)
  qv_scale(1) = str2num(a.as);
  qv_scale(2) = str2num(a.ds);
elseif ~isempty(a.as) && isempty(a.ds)
  qv_scale = [str2num(a.as) str2num(a.as)];
elseif isempty(a.as) && ~isempty(a.ds)
  qv_scale = [str2num(a.ds) str2num(a.ds)];
else
  qv_scale = [];
end

if ~isempty(a.aLabels) && ~isempty(a.dLabels)
  aLabels = a.aLabels;
  dLabels = a.dLabels;
elseif ~isempty(a.aLabels) && isempty(a.dLabels)
  aLabels = a.aLabels;
  dLabels = aLabels;
elseif isempty(a.aLabels) && ~isempty(a.dLabels)
  aLabels = a.dLabels;
  dLabels = aLabels;
else
  aLabels = []; dLabels = [];
end

if ~isempty(a.aTopLabels) && ~isempty(a.dTopLabels)
  aTopLabels = a.aTopLabels;
  dTopLabels = a.dTopLabels;
elseif ~isempty(a.aTopLabels) && isempty(a.dTopLabels)
  aTopLabels = a.aTopLabels;
  dTopLabels = aTopLabels;
elseif isempty(a.aTopLabels) && ~isempty(a.dTopLabels)
  aTopLabels = a.dTopLabels;
  dTopLabels = aTopLabels;
else
  aTopLabels = []; dTopLabels = [];
end

fname = [ num2str(size(xCL21.dat,2)) a.ext]

gistic_plots(base_dir,fname,xCL21,q,ads,regs,cyto,a.gn,qv_scale,aLabels,dLabels,aTopLabels,dTopLabels,a.qt,a.ull,a.wp,a.ns,a.abb)

param_file = [base_dir 'parameter_' method_st, '.txt'];
disp(['Writing parameter file to:' param_file]);
local_copy=1;
param_struct = struct('input_CL21_name', a.uiCL21, 'input_GISTICcore_filename', ...
                      a.GISTIC, 'cyto_file', a.uiCyto, 'parameter_file', ...
                      a.p ,'output_dir', base_dir, 'ext', a.ext ,'qthreshold', ...
                      a.qt,'ull', a.ull,'wp', a.wp,'ns', a.ns,'abb', a.abb, ...
                      'as', a.as,'ds', a.ds,'aLabels',a.aLabels, ...
                      'dLabels',a.dLabels, 'aTopLabels',a.aTopLabels, ...
                      'dTopLabels',a.dTopLabels,'genenames', a.gn);

gp_write_params(method_st,param_struct,params,param_file,local_copy);
