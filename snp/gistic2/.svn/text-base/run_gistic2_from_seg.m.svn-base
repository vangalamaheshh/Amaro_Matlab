function [D focal_regs params] = run_gistic2_from_seg(base_dir,seg_file,markers_file,ref_gene_file,params)
%RUN_GISTIC2_FROM_SEG run GISTIC (version 2.0) from a segmented data file
%
%  [D focal_regs params] = run_gistic2_from_seg(base_dir,seg_file,...
%                                  markers_file,ref_gene_file,params)
%
%  RUN_GISTIC2_FROM_SEG is a top level GISTIC pipeline function that
%  performs a GISTIC analysis on segmented copy-number data input stored in
%  a tab-delimited text file. The data are processed into a D structure
%  which is then passed to RUN_GISTIC20.
%
%  BASE_DIR is the output directory where GISTIC analyses are written.
%
%  SEG_FILE names a segmented data input file. Briefly, this is a
%  tab-delimited text file with a row for each segment and columns (1)
%  sample, (2) chromosome, (3) start base, (4) end base, (5) number of
%  markers, and (6) log2 ratio copy number. The file has a heaader row.
%
%  MARKERS_FILE names a file locating data markers. Briefly, this is a
%  tab-delimited text file without headers where each row specifies the
%  (1) marker name, (1) chromosome, and (3) position on the chromosome.
%
%  REF_GENE_FILE is the path to to a matlab-formatted reference gene file.
%
%  PARAMS is a structure containing optional parameters:
%
%   ** PARAMS.cnv_file specifies a file containing the genomic locations of
%    germline copy number variants which should be excluded from the analysis.
%   ** PARAMS.join_segment_size specifies the minimum size for a segment:
%    segments smaller than this value are combined to their largest
%    neighboring segment - the copy number of the resulting segment is a
%    weighted average. The default value is 8.
%   ** PARAMS.remove_X is a  boolean that removes chromosomes numbered
%    above 22 from the analysis. The default is 0
%   ** use_segarray is a boolean controlling memory performance: if 1, then
%    data is stored in a compressed format that takes more time to access.
%    The default is 0.
%   ** PARAMS.run_broad_analysis, if set, runs a broad analysis
%    (default cleared)
%   ** PARAMS.write_gene_files, if set, enables the output of
%    broad, focal, and combined broad/focal data by genes.
%   ** PARAMS.gene_collapse_method - string designating how to collapse marker
%    level data to gene level data for the gene files. 'mean' specifies to
%    use the average copy level of the markers (default); 'median' specifies
%    to  use the median marker value; 'extreme' means use the value furthest
%    from zero.
%   ** PARAMS.ext is a string specifying  a penultimate extension
%    added to saved file names to distinguish between runs. It
%    should start with, but not end with, a period. 
%   ** PARAMS.do_gene_gistic, if set, does gene_gistic for the
%    deletion analysis (default cleared)
%   ** PARAMS.t_amp and PARAMS.t_del are log2 copy number values
%    used to filter amplifications and deletions, and are also
%    passed to run_focal_gistic and gistic_broad_analysis.
%   ** PARAMS.broad_len_cutoff is the length cutoff used to
%    distinguish between broad and focal events. 
%   ** PARAMS.conf_level specifies the confidence levels of a peak
%    containing a driver used for setting the peak boundaries. If
%    passed a vector of values, GISTIC2 will report multiple peaks
%    for each region, one for each confidence level.  The default
%    confidence value is .75 (e.g. 75% confidence).
%   ** PARAMS.alpha specifies the exponential [amplification,deletion]
%    factors used to model the frequency dependence on event amplitude. The
%    default values are [2.5145 2.1653].
%   ** PARAMS.cap is a limit that is imposed on the event amplitudes
%    (+/- log ratio units). The default value is 1.5.
%   ** PARAMS.res specifies the resolution used to construct the empirical
%    frequency distributions used to calculate the background frequency.
%    Smaller values are more accurate, but result in longer computation
%    times. The default value is 0.05.
%   ** PARAMS.qv_thresh specifies the minimum q-value that a peak must have
%    to be considered significant. The default is .25.
%   ** PARAMS.array_list_file names a tab-delimited file with headers whose
%    first column contains the list of arrays on which to run the analysis.
%    If empty, analyzes all samples in D. 
%   ** PARAMS.do_arbitration is a boolean value: if 1 use arbitration
%    algorithm to distribute score between overlapping peaks (the default);
%    if 0 use greedy algorithm to give overlapping scores to the biggest peak.
%   ** PARAMS.arm_peeloff peel off just segments (0, the default) or segments
%    along with other events on the same arm for a given sample (1)
%    if 0 use greedy algorithm to give overlapping scores to the biggest peak.
%   ** PARAMS.peak_types specifies the type of wide peaks to calculate: 
%    'robust' for robust algorithm, 'loo' for leave-one-out algorithm.
%   ** PARAMS.use_two_sided if true, generate 2D quadrant figure for GISTIC
%    broad analysis. The default is false.
%   ** PARAMS.save_data_files is a boolean which controls the saving of
%    matlab data files to persistant storage. The default is 1.
%   ** PARAMS.conserve_disk_space is a boolean value which further
%    controls the data saving behavior: if 1, than the raw D structure,
%    focal segments and marker scores are not saved to persistant storage.
%   ** PARAMS.ziggs - is a structure containing ziggurat analysis parameters:
%      -- PARAMS.ziggs.max_segs_per_sample is the maximum number of
%        segments a sample can have before it is filtered out of the
%        analysis.
%   ** PARAMS.gene_collapse_method - string designating how to collapse marker
%    level data to gene level data. 'mean' specifies to use the average copy
%    level of the markers (default); 'median' means use the median marker value; 
%    'extreme' means use the value furthest from zero.


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  
  %% Set defaults  
  
  if ~exist('params','var') || isempty(params)
    params = struct();
  end
  
  % use standard GISTIC2 defaults for parameters that are not defined
  params = gistic2_param_defaults(params);
  % impose default values that are specific to using segmented data input
  params = impose_default_value(params,'cnv_file','');
  params = impose_default_value(params,'join_segment_size',8);
  params = impose_default_value(params,'remove_X',0);
  params = impose_default_value(params,'save_seg_data',1);
  params = impose_default_value(params,'islog',true);
  
  %% Check required inputs

  % base directory
  verbose('Checking inputs...',20);
  if ~exist('base_dir','var') || isempty(base_dir)
    error('Must supply base_dir!');
  else
    base_dir = add_slash_if_needed(base_dir);  
    if ~exist(base_dir,'dir')
      mkdir(base_dir)
    end
  end
  
  % segmentation file
  if ~exist('seg_file','var') || isempty(seg_file)
    error('Must supply seg file');
  else
    if ~exist(seg_file,'file')
      error('Must supply valid seg file!')
    end
  end
  
  % markers file
  if ~exist('markers_file','var') || isempty(markers_file)
    error('Must supply markers file');
  else
    if ~exist(markers_file,'file')
      error('Must supply valid markers file!')
    end
  end
  
  % reference genome file
  if ~exist('ref_gene_file','var') || isempty(ref_gene_file)
    error('Must supply ref gene file')
  else
    if ~exist(ref_gene_file,'file')
      error('Must supply valid ref gene file!')
    end
  end
  
  %% save inputs for future reproduction
  version = gistic_version; % LH string, RH function returns string
  save([base_dir 'gistic_inputs.mat']);

  %% Prepare D file
  verbose('Making D struct',20);
  options = struct('use_segarray',params.use_segarray,...
                   'islog',params.islog);
  CL = make_D_from_seg(seg_file,markers_file,options);
  verbosedisp(CL,30)
  
  %% Remove CNVs
  verbose('Removing CNVs',20)
  if ~isempty(params.cnv_file)
    [CL,nremoved] =remove_cnv(CL,params.cnv_file);
    fprintf(1,'Removed %d CNVs from raw data\r\n',nremoved);
  end
  verbosedisp(CL,30)
  
  % check that there are data left to process
  if size(CL.dat,1) == 0
      throw(MException('snp:run_gistic2_from_seg:AllDataRemoved',...
                       'All input data were removed after germline CNV processing.'));
  end
  
  %% match to array list file before removing NaNs 
  if ~isempty(params.array_list_file)
    verbose('Sub-selecting samples...',20)
    AL = read_array_list_file(params.array_list_file);
    use_arrays = {AL.array};
    [Mt,m1,m2]=intersect(use_arrays,CL.sdesc);
    if ~isunique(m1) || ~isunique(m2)
      error('either array list or samples are not unique');
    end
    if length(m1)==length(use_arrays)
      verbose(['Matched all ' num2str(length(m1)) ' samples in array list file'],10);
    else
      verbose(['Matched ' num2str(length(m1)) ' arrays out of ' num2str(length(use_arrays)) ' in array list file'], ...
              10);
    end
    CL=reorder_D_cols(CL,m2);
    CL=rmfield_if_exists(CL,'orig');
  else
    use_arrays = {};
  end
  verbosedisp(CL,30)
  verbose(['Matrix size ' num2str(size(CL.dat)) ],20);
  % check that there are data left to process
  if size(CL.dat,1) == 0
      throw(MException('snp:run_gistic2_from_seg:AllDataRemoved',...
                       'All input data were removed after array list processing.'));
  end

  %% Remove NaN probes
  verbose('Removing NaN probes...',20)
  nan_idx=find(any(isnan(CL.dat),2));
  if ~isempty(nan_idx)
    verbose(['Removing ' num2str(length(nan_idx)) ' markers with NaNs'],20); 
    CL=reorder_D_rows(CL,setdiff(1:size(CL.dat,1),nan_idx));
  else
    verbose('No markers with NaNs... ',20);
  end  
  verbose(['Matrix size ' num2str(size(CL.dat)) ],20);
  % check that there are data left to process
  if size(CL.dat,1) == 0
      throw(MException('snp:run_gistic2_from_seg:AllDataRemoved',...
                       'All input data were removed after NaN processing.'));
  end
  
  %% Check if data is log, if not, take log
  if ~isfield(params,'islog')
    verbose('Checking if data is log2',20);
    s = nanmean(CL.dat,1);
    s = mean(s);
    if abs(s) >= .5
      verbose('Computing Log2 of Data',20);
      CL.dat = log2(CL.dat)-1;
    end
    CL.islog = 1;
  end
  %% remove X,Y chromosome
  if params.remove_X
    verbose('Removing X chromosome...',20);
    CL=reorder_D_rows(CL,CL.chrn<=22);
  end
  CL=rmfield_if_exists(CL,{'cbs','cbs_rl'});
  
  %% Join small segments
  CL.cbs=CL.dat;
  verbose('Merging small segments...',20);
  CL=smooth_cbs(CL,params.join_segment_size);
  CL.dat=CL.cbs;
  CL=rmfield_if_exists(CL,{'cbs','cbs_rl'});
  verbosedisp(CL,30)
  
  %% subtract median
  verbose('Median centering data...',20);
  CL.medians=median(CL.dat((CL.chrn<=22),:),1);
  CL.dat=CL.dat-repmat(CL.medians,size(CL.dat,1),1);
  
  %% save segmented data
  if params.save_seg_data && ~params.genepattern
    verbose('Saving segmented data file...',20);
    seg_data_file = [base_dir params.fname 'segmented_data' params.ext '.mat'];
    save_D(seg_data_file,CL,'-v7.3');  
  end
  
  %% Now Run GISTIC 2.0
  [D focal_regs params] = run_gistic20(base_dir,CL,ref_gene_file,params);
  
 
  
