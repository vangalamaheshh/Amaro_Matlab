ptable=load_struct('/Volumes/xchip_cga_home/amaro/TumorInNormal/Pancan_NT_paper/pairtable.tsv');
for i=1:slength(ptable)
    ptable.individual_id{i,1}=ptable.pair_id{i}(1:12);
end
for i=1:slength(ptable)
    ptable.sample_type{i,1}=ptable.pair_id{i}(end-1:end);
end
[m l]=count(ptable.individual_id);
[i n]=ismember(ptable.individual_id,l);
ptable.pair_count(i)=m(n);

analysis_table=reorder_struct(ptable,~ismember(ptable.call_stats_capture,'')&ptable.pair_count>1);

[m l]=count(analysis_table.individual_id);
[i n]=ismember(analysis_table.individual_id,l);
analysis_table.pair_count(i)=m(n);
analysis_table=reorder_struct(analysis_table,analysis_table.pair_count>1);