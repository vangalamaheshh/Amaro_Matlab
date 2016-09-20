function run_cbs_batch(CN,ch,should_log,should_smooth,window_size,samples_idx)

warning(['This function has been deprecated. Please use run_cbs ' ...
         'instead. NS 2009.02.09'] )

if ~exist('samples_idx','var')
  samples_idx = [1:getsize(CN,'dat',2)];
end
  
for i=1:length(samples_idx)
    run_cbs(['Sample' sprintf('%03d',samples_idx(i))],reorder_D_cols(CN,samples_idx(i),'allmem'),ch,should_log,should_smooth,window_size);
end