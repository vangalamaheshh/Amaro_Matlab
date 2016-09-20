for i=1:slength(individual_table)
if sum(ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'04')&ismember(sample_table.sample_type,'Tumor'))>0
    individual_table.tp4{i,1}=sample_table.sample_id{ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'04')&ismember(sample_table.sample_type,'Tumor')};
else
    individual_table.tp4{i,1}='';
end
if sum(ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'03')&ismember(sample_table.sample_type,'Tumor'))>0
individual_table.tp3{i,1}=sample_table.sample_id{ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'03')&ismember(sample_table.sample_type,'Tumor')};
else
    individual_table.tp3{i,1}='';
end
if sum(ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'02')&ismember(sample_table.sample_type,'Tumor'))>0
individual_table.tp2{i,1}=sample_table.sample_id{ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'02')&ismember(sample_table.sample_type,'Tumor')};
else
    individual_table.tp2{i,1}='';
end
if sum(ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'01')&ismember(sample_table.sample_type,'Tumor'))>0
individual_table.tp1{i,1}=sample_table.sample_id{ismember(sample_table.individual,individual_table.individual_id{i})&ismember(sample_table.time_point,'01')&ismember(sample_table.sample_type,'Tumor')};
else
    individual_table.tp1{i,1}='';
end
end