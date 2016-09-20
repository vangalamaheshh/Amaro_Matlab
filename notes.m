for j=1:slength(homs)                          
h=load_struct(homs.col1{j});
for i=1:slength(CallStats)-1                   
c=load_struct(CallStats.call_stats_capture{i});
x(j,i)=sum(ismember(c.position,h.position));
calls{i}=CallStats.pair_id{i};
end
homs_names{j}=h.Sample{1};
j
end