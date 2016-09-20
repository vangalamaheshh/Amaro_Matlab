function gp_gistic2_from_seg(varargin)
%GP_GISTIC2_FROM_SEG gene pattern module for GISTIC 2.0
%
%   gp_gistic2_from_seg -b base_dir -seg segfile  -mk markersfile 
%   -refgene refgenefile [-alf array_list_file] [-cnv cnv_file] 
%   [-ta t_amp] [-td t_del] [-js join_segment_size] [-qvt qv_thresh]
%   [-rx remove_X] [-v verbose_level] [-ext extension] [-cap cap_val]
%   [-broad run_broad_analysis] [-brlen broad_length_cutoff]
%   [-maxseg max_sample_segs] [-res res] [-conf conf_level] 
%   [-genegistic do_gene_gistic] [-smalldisk save_disk_space]
%   [-smallmem use_segarray] [-savegene write_gene_files]
%   [-arb do_arbitration] [-twosides use_two_sided] [-peaktype peak_types]
%   [-saveseg save_seg_data] [-savedata write_data_files]
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


%  COMPILE INSTRUCTIONS: (as of 2011-03-22)
%       1) dotkit environment
%           a) use Matlab-2010b  (MATLAB v7.14, compiler v4.14)
%           b) use GCC-4.3  (minimum C/C++ compiler for this Matlab version)
%           c) use Matlab_MCR (locates the Matlab Component Runtime: a
%              caveat is that this will point to the wrong runtime once
%              R2010b is no longer current).
%
%       2) Need to run 'mbuild -setup' at least once for the compiler to
%          work
%
%       3) change to the output directory (where the executable will be emitted)
%          and type:
%              mcc -v -m -w enable gp_gistic2_from_seg
%
%---
% $Id$
% $Date: 2011-10-05 18:08:08 -0400 (Wed, 05 Oct 2011) $
% $LastChangedBy: schum $
% $Rev$

% create documentation-only dependencies
Qs;
% OpenGL hardware rendering is sometimes problematic
opengl software

if isempty(varargin)
    usage();
    return
end

a = handle_args({'b','seg','mk','refgene','alf','cnv','ta','td','js',...
                'ext','qvt','rx','v','broad','brlen','maxseg','res',...
                'conf','cap','genegistic','smalldisk','savedata','smallmem',...
                'savegene','twoside','arb','peaktype','saveseg','genepattern'},varargin);

% validate arguments
try
    [base_dir segfile markersfile refgenefile] = validate_mandatory_args(a);
    optional_params = validate_optional_args(a);
catch me
    if strmatch('snp:',me.identifier)
        % intercept and report user input errors
        disp(' ');
        disp('GISTIC 2.0 input error:');
        disp(me.message);
        disp(' ');
        usage();
        return
    else
        % report errors we did not deliberately generate for debugging
        rethrow(me);
    end
end


%% run GISTIC 2.0 from segmented data

disp(['GISTIC version ' gistic_version]);

% display parameters
verbose('required parameters:',10)
verbose(['     base_dir: ' base_dir],10);
verbose(['      segfile: ' segfile],10);
verbose(['  markersfile: ' markersfile],10);
verbose(['  refgenefile: ' refgenefile],10);
verbosedisp(optional_params,10);

% run gistic in try-catch block
try
    run_gistic2_from_seg(base_dir,segfile,markersfile,refgenefile,optional_params);
catch me
    if strmatch('snp:',me.identifier)
        % intercept and report user input errors
        disp(' ');
        disp('GISTIC 2.0 input error:');
        disp(me.message);
        return
    else
        % report errors we did not deliberately generate for debugging
        rethrow(me);
    end
end
        
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% usage message
function usage
    disp('Usage: gp_gistic2_from_seg -b base_dir -seg segmentation_file -mk markers_file');
    disp('-refgene ref_gene_file [-alf array_list_file(def:empty)]');
    disp('[-cnv cnv_file] [-ta amplifications_threshold(def=.1)] [-td deletions_threshold(def=.1)]');
    disp('[-js join_segment_size(def=8)] [-ext extension] [-qvt qv_thresh(def=0.25)]');
    disp('[-rx remove_x(def=1)] [-v verbosity_level(def=0)] [-cap cap_val(def=1.5]]');
    disp('[-broad run_broad_analysis(def=0)] [-brlen broad_length_cutoff(def=0.98)]');
    disp('[-maxseg max_sample_segs(def=2500)] [-res res(def=0.05)] [-conf conf_level(def=0.75)]');
    disp('[-genegistic do_gene_gistic(def=0)] [-smalldisk save_disk_space(def=0)]');
    disp('[-smallmem use_segarray(def=1)] [-savegene write_gene_files(def=0)]');
    disp('[-arb do_arbitration(def=1)] [-twosides use_two_sided(def=0)] [-peaktype peak_types(def=robust)]');
    disp('[-saveseg save_seg_data(def=1)] [-savedata write_data_files(def=1)]');

%% validate mandatory arguments
%
function [base_dir segfile markersfile refgenefile] = validate_mandatory_args(a)
% validate base output directory
base_dir = a.b;
if isempty(base_dir)
    throw(MException('snp:badarg','Must supply a base directory.'))
end
if ~exist(base_dir,'dir')
    throw(MException('snp:badarg','Base directory ''%s'' does not exist',a.b));
end

% validate segmentation file
segfile = a.seg;
if isempty(segfile)
    throw(MException('snp:badarg','Must supply a segmentation file.'));
end
if ~exist(segfile,'file')
    throw(MException('snp:badarg','Segmentation file ''%s'' does not exist',segfile));
end

% validate markers file
markersfile = a.mk;
if isempty(markersfile)
    throw(MException('snp:badarg','Must supply a markers file.'));
end
if ~exist(markersfile,'file')
    throw(MException('snp:badarg','Markers file ''%s'' does not exist',markersfile));
end

% validate refgene file
refgenefile = a.refgene;
if isempty(refgenefile)
    throw(MException('snp:badarg','Must supply a refgene file.'));
end
if ~exist(refgenefile,'file')
    throw(MException('snp:badarg','Reference gene file ''%s'' does not exist',refgenefile));
end


%% validate optional parameters
% evaluate a struct of user command line settings for known parameters
% return a struct of gistic optional params
function params = validate_optional_args(a)

%% gistic 1.0 optional arguments
params = struct;

if ~isempty(a.alf)
    params.array_list_file = a.alf;
else
    params.array_list_file = [];
end

if ~isempty(a.cnv)
    params.cnv_file = a.cnv;
else
    params.cnv_file = [];
end

% amplification and deletion thresholds
params.t_amp = numeric_arg(a,'ta',0.1,[0 Inf]);
params.t_del = numeric_arg(a,'td',0.1,[0 Inf]);

% join segment size
params.join_segment_size = numeric_arg(a,'js',8,[1 Inf]);

% output file extension
if isempty(a.ext)
    params.ext=[];
else
    params.ext=a.ext;
end

% q-value threshold
params.qv_thresh = numeric_arg(a,'qvt',0.25,[0 1]);
% remove X
params.remove_X = numeric_arg(a,'rx',1,[0 1]);

if ~isempty(a.v)
    set_verbose_level(str2double(a.v));
end

%% gistic 2.0 optional parameters

% broad analysis flag
params.run_broad_analysis = numeric_arg(a,'broad',0,[0 1]);
% broad/focal length cutoff
params.broad_len_cutoff = numeric_arg(a,'brlen',0.98,[0 2]);
% maximum segments per sample
maxseg = numeric_arg(a,'maxseg',2500,[0 Inf]);
params.ziggs = struct('max_segs_per_sample',maxseg);

params.res = numeric_arg(a,'res',0.05,[0 1]);
params.conf_level = numeric_arg(a,'conf',0.75,[0 1]);
params.cap = numeric_arg(a,'cap',1.5,[0 Inf]);
params.do_gene_gistic = numeric_arg(a,'genegistic',0,[0 1]);
params.conserve_disk_space = numeric_arg(a,'smalldisk',0,[0 1]);
params.save_data_files = numeric_arg(a,'savedata',1,[0 1]);
params.use_segarray = numeric_arg(a,'smallmem',1,[0 1]); % note: overrides source default
params.write_gene_files = numeric_arg(a,'savegene',0,[0 1]);
params.use_two_sided = numeric_arg(a,'twoside',0,[0 1]);
params.do_arbitration = numeric_arg(a,'arb',1,[0 1]);
params.save_seg_data = numeric_arg(a,'saveseg',1,[0 1]);
% validate peaktype argument
if ~isempty(a.peaktype)
    if any(strcmp(a.peaktype,{'robust','loo'}))
        params.peak_types = a.peaktype;
    else
        throw(MException('snp:badarg','-peaktype argument must be ''robust'' or ''loo''.'));
    end
else
    params.peak_types = 'robust';
end
% make sure peak_types is a cell array
if ~iscell(params.peak_types)
  params.peak_types = {params.peak_types};
end
params.genepattern = numeric_arg(a,'genepattern',1,[0 1]);

%% numeric parameter handler
%
function val = numeric_arg(args,param,default,range)
    if isempty(args.(param))
        val = default;
    else
        val = str2double(args.(param));
        if isnan(val)
            throw(MException('snp:badarg:nonnumeric','Non-numerical value ''%s'' for -%s',...
                     args.(param),param));
        end
    end
    % optional range check
    if exist('range','var') && ~isempty(range)
        if val < range(1) || val > range(2)
            throw(MException('snp:badarg:outofrange',...
                    ['Value of ' num2str(val) ' for -' param ' is outside of allowed range ' ...
                    num2str(range(1)) ' - ' num2str(range(2)) '.']));
        end
    end
