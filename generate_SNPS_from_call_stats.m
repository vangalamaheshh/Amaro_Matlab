function generate_SNPS_from_call_stats(germline_snps,pair_id)
call=load_struct(germline_snps); % PoN annotated call_stats file required
disp('loaded call stats file preprocessing file')
call.contig=chromosome2num_legacy(call.contig);
call.position=str2double(call.position);
call.t_alt_count=str2double(call.t_alt_count);
call.t_ref_count=str2double(call.t_ref_count);
call.n_alt_count=str2double(call.n_alt_count);
call.n_ref_count=str2double(call.n_ref_count);
call.PoN_Germline=str2double(call.PoN_Germline);
call.alt_count_greater10_af_greater_20percent=str2double(call.alt_count_greater10_af_greater_20percent);
call.normal_f=str2double(call.normal_f);
    call.tumor_f=str2double(call.tumor_f);




if isfield(call,'total_reads')
call.total_reads=str2double(call.total_reads);
else
    call.total_reads=call.t_alt_count+call.t_ref_count+call.n_alt_count+call.n_ref_count;
end
call.observed_in_normals_count=str2double(call.observed_in_normals_count);


af_range=[0:.01:1];
call.PoN_Artifact=str2double(call.PoN_Artifact);

xCoV=call.total_reads>10;
xPoN=(call.alt_count_greater10_af_greater_20percent>.05&call.observed_in_normals_count>30);
xART=call.PoN_Artifact<.1;
xAF=call.normal_f<.7&call.normal_f>.3;
call=reorder_struct(call,xPoN&xAF&xART&xCoV);



call=reorder_struct(call,~isnan(call.contig));
call=reorder_struct(call,~(call.contig==23|call.contig==24));
save_struct(call,strcat(pair_id,'.hetsites.txt'));
quit
end