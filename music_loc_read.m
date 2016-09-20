d = dataset('file', '~/Downloads/artists.csv', 'ReadObsNames', false, 'ReadVarNames', true, 'Delimiter', ',');
%artist_location=cell(size(d,1),1);


parfor k = 1 : size(d,1)
    k
    error=[];
     try x=urlread(['https://musicbrainz.org/artist/', d.artist{k}]);
     catch error
     end
     if isempty(error)
    loc_string=regexp(x,'<dd class="area">[/\\,[0-9],<,>,",a-z,\s,=,[A-Z],-]+</dd>','match');
    if ~isempty(loc_string)
    loc_string=regexp(loc_string,'<bdi>[[,],A-Z,a-z,\s]+</bdi>','match');
    loc_string=regexprep(loc_string{1},'<bdi>','');
    artist_location(k,1)=regexprep(loc_string(end),'</bdi>','');
    end
     end
end