

cd ~/Projects/Colorectal/
calls_tsv=load_struct('raw_normal_call.tsv');
bed=load_struct('CRC_starter_genes.txt');
bed.start=str2double(bed.start);
bed.end=str2double(bed.end);

call=load_struct(calls_tsv.call_stats_capture_germline{1});
call.position=str2double(call.position);
fishing=fish_for_genes(call,bed);

for i=1:slength(calls_tsv)
    if ~isempty(calls_tsv.call_stats_capture_germline{i})
    call=load_struct(calls_tsv.call_stats_capture_germline{i});
    call.position=str2double(call.position);
    call=reorder_struct(call,ismember(call.contig,{'3';'5';'12';'17'}));
    fishing_in=fish_for_genes(call,bed);
    fishing=mergeStruct(fishing,fishing_in);
    if mod(i,20)==0
        disp(sprintf('Searching calls file %d',i))
        count(fishing.contig)
    end
    end
end

