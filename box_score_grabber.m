function [ box_scores handle_web teamnames] = box_score_grabber( URL, numberofgames, i)
%This function takes a URL and grabs the boxscores found there
%For now this assumes basketball reference format as of 8/14/2012


[~, handle_web url]=web(URL); pause(10);

getTableFromWeb();
pause(2);
box_scores{1}=getTableFromWeb(numberofgames+4); 
box_scores{2}=getTableFromWeb(numberofgames+5); 
box_scores{3}=getTableFromWeb(numberofgames+6);
box_scores{4}=getTableFromWeb(numberofgames+7);


teamnames{1}=getTableFromWeb(i);

end

