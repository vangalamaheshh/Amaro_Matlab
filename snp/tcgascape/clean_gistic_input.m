function D = clean_gistic_input(D,cnv_file,remove_XYM,join_segment_size)

% remove areas of normal copy number variation
if ~isempty(cnv_file)
    D = remove_cnv(D,cnv_file);
end

% Remove NaN probes
verbose('Removing NaN probes...',20)
nan_idx=find(any(isnan(D.dat),2));
if ~isempty(nan_idx)
    verbose(['Removing ' num2str(length(nan_idx)) ' markers with NaNs'],20); 
    D=reorder_D_rows(D,setdiff(1:size(D.dat,1),nan_idx));
else
    verbose('No markers with NaNs... ',20);
end  
verbose(['Matrix size ' num2str(size(D.dat)) ],20);
  
% remove X,Y chromosome
if remove_XYM
    verbose('Removing X chromosome...',20);
    D=reorder_D_rows(D,D.chrn<=22);
end

% join small segments 
verbose('Merging small segments...',20);
D = rmfield_if_exists(D,{'cbs','cbs_rl'});
D.cbs=D.dat;
D = smooth_cbs(D,join_segment_size);
D.dat=D.cbs;
D = rmfield_if_exists(D,{'cbs','cbs_rl'});

% subtract median
verbose('Median centering data...',20);
D.medians = median(D.dat,1);
D.dat = D.dat-repmat(D.medians,size(D.dat,1),1);
