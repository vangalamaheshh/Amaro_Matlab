function HetSitesFilterOnFrequency(calls,freq)
addpath ../../../../matlab/
addpath ../../../../matlab/mike/
addpath ../../../../matlab/seq/


load(freq);
calls=load_table(calls);

calls=rmfield(calls,'header');
calls=rmfield(calls,'headline');

calls=reorder_struct(calls,ismember(calls.position,pop_numbers.position));
save_struct(calls,'call.tmp');


end
