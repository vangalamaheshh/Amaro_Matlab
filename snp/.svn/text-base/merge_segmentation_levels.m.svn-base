function C1=merge_segmentation_levels(C,merge_type)

if ischar(merge_type)
  tmp.method=merge_type;
  merge_type=tmp;
end

C=remove_X(C);
% C=sub_median_from_raw(C,0); % cbs_rl is recalculated from median fixed values
C=calc_segment_std_err(C,0);
%% FIXME: median of raw does NOT fit cbs_rl
%% at least the cbs_rl and dat represent the same data

switch merge_type.method 
 case 'histogram'
%  disp([ min(C.dat(:)) max(C.dat(:)) ]);
  C1=hist_of_smoothed(C,(min(C.dat(:))-3*0.03):0.001:(max(C.dat(:))+3*0.03),0.03);
 case 'cluster'
  if isfield(merge_type,'lsf_params')
    C1=merge_levels(C,merge_type.th,merge_type.lsf_params);% -log(1e-10)
  else
    C1=merge_levels(C,merge_type.th);% -log(1e-10)
  end    
end

