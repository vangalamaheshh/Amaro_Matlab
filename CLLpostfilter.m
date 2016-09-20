
call_stats_files_est=load_struct('~/Projects/CLL_Rush/ABSOLUTE/FilteredMafs/22MAFs/Mafs');
call_stats_files_est.tumor_in_normal=str2double(call_stats_files_est.tumor_in_normal);
for i=1:slength(call_stats_files_est)
calls=load_struct(call_stats_files_est.call_stats_capture_no{i});
calls.t_alt_count=str2double(calls.t_alt_count);
calls.t_ref_count=str2double(calls.t_ref_count);
calls.n_alt_count=str2double(calls.n_alt_count);
calls.n_ref_count=str2double(calls.n_ref_count);

for f=1:slength(calls)
calls.logOdds(f,1)=normal_contamination_estimate_fitting(calls.t_alt_count(f),calls.t_ref_count(f)...
    ,calls.n_alt_count(f),calls.n_ref_count(f),call_stats_files_est.tumor_in_normal(i),.95);
end

ratio_filter=ismember(calls.judgement,'KEEP')&calls.logOdds<7;  
save_struct(calls,sprintf('~/Projects/CLL_Rush/ABSOLUTE/FilteredMafs/22MAFs/%s.maf',call_stats_files_est.pair_id{i}));
end