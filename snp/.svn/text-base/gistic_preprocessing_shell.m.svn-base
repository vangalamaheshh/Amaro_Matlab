function [D,repnames,bad_tumor_names,arrays_for_core,celllines] = gistic_preprocessing(base_dir,...
    sample_info_file,array_list_file,Dscell,...
    snp_skip,show_hist,save_raw,norm_collapse_method,use_paired,...
    n_closest_n,perform_batch_correction,batch_effect_min_sz,batch_effect_bonf_pv_thresh,...
    batch_effect_abs_pv,output_dir_extension,hist_qc_normals,hist_qc_tumors,...
    conserve_memory,memlimit,use_all_cores,tan_ceil_val,scaleoper,startfromraw,genepattern,cnvfile,genomebuild)  %#ok<INUSL,INUSD>
%GISTIC_PREPROCESSING preprocess array data.  (Batch effect correction,
%normalization, and hist_qc.)
%
%   DOUT = GISTIC_PREPROCESSING(BASE_DIR, SAMPLE_INFO_FILE,
%   ARRAY_LIST_FILE, DFILESCELL, SNP_SKIP, SHOW_HIST, SAVE_RAW,
%   NORM_COLLAPSE_METHOD, USE_PAIRED, N_CLOSEST_N,
%   PERFORM_BATCH_CORRECTION, BATCH_EFFECT_MIN_SZ,
%   BATCH_EFFECT_BONF_PV_THRESH, BATCH_EFFECT_ABS_PV, OUTPUT_DIR_EXTENSION,
%   HIST_QC_NORMALS, HIST_QC_TUMORS, CONSERVE_MEMORY,MEMLIMIT,use_all_cores,tan_ceil_val,
%   scaleoper,startfromraw,genepattern,cnvfile,genomebuild)
%
%
%   Description of Inputs:
%
%   base_dir: (required) base directory for data processing; (string)
%   sample_info_file: (required) sample info file (in base dir) (string)
%   array_list_file: name of array list file (can be cell array of strings(?)), default empty (cell array of strings)
%   dscell: cell array listing filenames of D struct files (from
%       snp_to_D); pathnames are relative to base directory  OR cell array
%       of structs
%   snp_skip: decimation frequency on snp data, default 1 (int)
%   show_hist: show histogram?, 0 = no (default); 1 = yes (logical)
%   save_raw: save raw data?, 0 = no (default); 1 = yes (logical)
%   norm_collapse_method: 'mean', 'median', or 'tangent' (string)
%   use_paired: use_paired_norms? 0 = no (default); 1 = yes 
%   n_closest_n: how many norms for n_closest_n? default = 5 (int)
%   perform_batch_correction: 0 = no (default); 1 = yes (logical)
%   batch_effect_min_sz: minimum size batch to correct, default = 5 (int)
%   batch_effect_bonf_pv_thresh: pre-Bonferroni correction threshold,
%        default=0.05 (double) **NOTE: BONFERRONI CORRECTION NEVER
%        HAPPENS**
%   batch_effect_abs_pv: p value threshold for no Bonferroni correction,
%        default = 0.001 (double)
%   output_dir_extension: output dir extension (default = <yymmdd>)
%        (string)
%   hist_qc_normals: use hist_qc to remove normals with copy number changes?
%        (0 = no (default); 1 = yes) (logical)
%   hist_qc_tumors: use hist_qc to remove tumors without copy number changes?
%        (0 = no; 1 = yes (default))  (logical)
%   conserve_memory=1: uses shorts and uint8 to save memory
%                  =2: converts dstructs from structs to datastructs so that
%   data and affy_calls is saved as HDF5 files on disk.
%        (0 = no (default))
%   memlimit: when implemented in a subfunction, the maximum number of
%   bytes of data from a field in datastruct object to load into memory at
%   a time (should use conserve_memory=2). Default: 300000000 (300MB).
%   use_all_cores: use all cores for reference
%   tan_ceil_val: default empty
%   scaleoper: 'mean' by default.  Also accepts function handles.
%   startfromraw: 1 = load D.raw.mat rather than individual plates (Dscell should then have location of D.raw.mat file)
%   genepattern: 1 if running as executable, 0 if running in matlab (won't
%   datastructure files)
%   cnvfile (optional) the name of the cnvfile 
%   genomebuild: {hg16,hg17,hg18}
%  

%REMOVED:
%   save_after_tumor_qc: no = 0 (default); yes = 1.  (Save data before
%   removing bad tumors and samples not intended for gistic core?)
%   maxsamps, maxsnps
%
%   Change on 01 Feb 2008:  Selecting tangent normalization does not automatically cause histqc to be run on normals first. 
%  
% ---
% $Id$
% $Date: 2007-12-03 15:57:57 -0500 (Mon, 03 Dec 2007) $
% $LastChangedBy: rameen $
% $Rev$



%% Check Inputs

varlist1 = {'base_dir','sample_info_file','array_list_file','Dscell','snp_skip',...
    'show_hist','save_raw','norm_collapse_method','use_paired',...
    'n_closest_n','perform_batch_correction','batch_effect_min_sz','batch_effect_bonf_pv_thresh',...
    'batch_effect_abs_pv','output_dir_extension','hist_qc_normals','hist_qc_tumors','conserve_memory',...
    'memlimit','use_all_cores','tan_ceil_val','scaleoper','startfromraw','genepattern','cnvfile','genomebuild'};

defaults = {'ERR','ERR','ERR','{}','1','0','0','''median''','0','5','1','5','.05',...
    '.001','''output''' ,'0','1','0','300000000','1','[]','''mean''','0','0','[]','''hg18'''};

required = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

for idx = 1:length(varlist1)
    if required(idx) && (~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')']))
        error('Required input %s undefined.',varlist1{idx})
    elseif ~exist(varlist1{idx},'var') || eval(['isempty(' varlist1{idx} ')'])
        eval([varlist1{idx} '=' defaults{idx} ';'])
    end
end

%% call gistic_prep function with options

% pack the optional parameters into an options struct
opts.snp_skip = snp_skip;
opts.show_hist = show_hist;
opts.save_raw = save_raw;
opts.norm_collapse_method = norm_collapse_method;
opts.use_paired = use_paired;
opts.n_closest_n = n_closest_n;
opts.perform_batch_correction = perform_batch_correction; %!
opts.batch_effect = struct('min_sz', batch_effect_min_sz,...
                           'bonf_pv_thresh',batch_effect_bonf_pv_thresh,...
                           'absolute_pv',batch_effect_abs_pv);

opts.output_dir_extension = output_dir_extension;
opts.hist_qc_normals = hist_qc_normals;
opts.hist_qc_tumors = hist_qc_tumors;
opts.conserve_memory = conserve_memory;
opts.memlimit = memlimit;
opts.use_all_cores = use_all_cores;
opts.tan_ceil_val = tan_ceil_val; %! needed? only use commented out
opts.scaleoper = scaleoper;
opts.startfromraw = startfromraw;
opts.genepattern = genepattern;
opts.cnvfile = cnvfile;
opts.cnv_byposition = 0;
opts = impose_default_value(opts,'genomebuild', 'hg18'); %! TODO should be file

% genomebuild has several hardwired path selections for hg16/17/18 
if strcmp('hg18',opts.genomebuild)
    rg_file = '/xchip/gistic/variables/hg18/hg18_20080402.mat';
elseif strcmp('hg17',genomebuild)
    rg_file = '/xchip/gistic/variables/cyto_rg_hg17_ucsc20070227.mat';
elseif strcmp('hg16',genomebuild)
    rg_file = '/xchip/gistic/variables/hg16_20070112.mat';
else
    rg_file = genomebuild;
end

% override some defaults
opts.log2_output = 0;

% call the inner function
[D,repnames,bad_tumor_names,arrays_for_core,celllines] = ...
    gistic_prep(base_dir,sample_info_file,array_list_file,Dscell,rg_file,opts);


