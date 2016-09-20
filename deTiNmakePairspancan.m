counter=1;
for i=1:length(indivs)
    T=find(ismember(table.short_letter_code,'TP')&ismember(table.individual_id,indivs{i}),1);
    
    NT=find(ismember(table.short_letter_code,'NT')&ismember(table.individual_id,indivs{i}),1);
    if ~isempty(T) && ~isempty(NT)
    pair_table.pair_id{counter,1}=strcat([table.sample_id{T},'_',table.sample_id{NT},'_',table.primary_disease{T},'_TP_NT']);
    pair_table.case_sample{counter,1}=table.sample_id{T};
    pair_table.control_sample{counter,1}=table.sample_id{NT};
    pair_table.individual_id{counter,1}=indivs{i};
    counter=counter+1;
if ismember('NB',table.short_letter_code(ismember(table.individual_id,indivs{i})))
    NB=find(ismember(table.short_letter_code,'NB')&ismember(table.individual_id,indivs{i}),1);
    pair_table.pair_id{counter,1}=strcat([table.sample_id{T},'_',table.sample_id{NB},'_',table.primary_disease{T},'_TP_NB']);
    pair_table.case_sample{counter,1}=table.sample_id{T};
    pair_table.control_sample{counter,1}=table.sample_id{NB};
    pair_table.individual_id{counter,1}=indivs{i};
    counter=counter+1;
end
    end
    
end
