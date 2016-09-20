function [result_out]=fish_for_genes(call,bed)

k=(call.position>bed.start(1)&call.position<bed.end(1)&ismember(call.contig,bed.chr{1}));
result_out=reorder_struct(call,k);
for i=2:slength(bed)
    k=(call.position>bed.start(i)&call.position<bed.end(i)&ismember(call.contig,bed.chr{i}));
    result=reorder_struct(call,k);
end

result_out=mergeStruct(result,result_out);
end