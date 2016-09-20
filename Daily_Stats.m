function [ stats_for_the_day teamnames ] = Daily_Stats(fake_date)
%This is the function that will run each day and grab the links for each
%box score on that day.
%It then takes the links for each day and runs box_score_grabber
%Will return with some var for the box_scores of each day

%fake_date=;
[month day year]=parse_date(fake_date);

url_games=sprintf('http://www.basketball-reference.com/boxscores/index.cgi?month=%s&day=%s&year=%s',month,day,year);

pageString=urlread(url_games);

[~, ~, tokenidx, matchstring]=regexp(pageString,'\d*[A-Z]*.html">Final');

if isempty(matchstring)
    disp('No Games Played on'); 
    month 
    day 
    year
   stats_for_the_day=[]; teamnames=[];

else
for i=1:size(matchstring,2)
    url_game=sprintf('http://www.basketball-reference.com/boxscores/%s',matchstring{i}(1:size(matchstring{i},2)-7));
    [stats_for_the_day{i} webhandle teamnames{i}]=box_score_grabber(url_game,size(tokenidx,2),i);
end

%close(webhandle)
clear mex

end
end
