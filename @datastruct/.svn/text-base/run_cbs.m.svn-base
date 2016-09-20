function C2 = run_cbs(prefix,C,ch,should_log,should_smooth,window_size,...
            rewrite_data_if_exists,write_header_for_output_file, localLSF)
%RUN_CBS run circular binary segmentation on a datastruct class object
%
%   C2 = run_cbs(PREFIX,C,CH,SHOULD_LOG,SHOULD_SMOOTH,WINDOW_SIZE,...
%            REWRITE_DATA_IF_EXISTS,WRITE_HEADER_FOR_OUTPUT_FILE, LOCALLSF)
%
% 
%

if exist('localLSF','var') && strcmp(localLSF,'local')
    warning('Local processing not supported for HDF5. Use LSF instead.');
end

if exist('should_smooth','var') && should_smooth
    warning('Smoothing not supported for HDF5.');
end

%%%%%% VIGNETTE
params = struct;
if ~exist('window_size','var') || ~isempty(window_size)
    params.should_window = 0;
else
    params.should_window = 1;
    params.window_size = window_size;
end

% separate path to directory from actual prefix
pathend = find(prefix==filesep,1,'last');
if isempty(pathend)
    pathend = -1;
    seg_dir = prefix(1:pathend);
else
    seg_dir = './';
end
prefix = prefix(pathend+1:end);
%run HDF5 segmentation function
run_cbs_batch_hdf5(C,[],seg_dir,params,[],~should_log,prefix(pathend+1:end));
C2 = []; % don't ask