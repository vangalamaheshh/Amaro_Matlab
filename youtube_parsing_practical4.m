
d = dataset('file', '~/Downloads/artists.csv', 'ReadObsNames', false, 'ReadVarNames', true, 'Delimiter', ',');
view_count=zeros(size(d,1),10);
parfor k = 1 : size(d,1)
    k
    art = regexprep(d.name{k}, '\s+', '+');
if isequal(art,'!!!')
art='bangbangbang';
end
    x=urlread(['https://www.youtube.com/results?search_query=', art]);
    match=regexp(x,'>[0-9,\,]+ views<','match');
    match=char2cell(match);
    for i=1:10
        view_count(k,i)=parse_match_string(match{1}{i});
    end
end


