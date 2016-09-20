%% Fantasy Basketball League Machine%%
%%% NERD BASKETBALL TIME
%% generates a spreadsheet of the boxscores of the day with player names and stats
%Teams=readtext('Teams.txt');
% Teams{1}={'Andre Iguodala'};
% player_loc={};

for day_n=1:size(fake_date,2)
keep('day_n', 'fake_date')



[Stats Team_Names] =Daily_Stats(char(fake_date{day_n}));
if isempty(Stats)
else

 T_Codes=cell(size(Team_Names,2),2);
 for n=1:size(Team_Names,2)
     for i=1:size(Team_Names{n})%should always be 1
         %for j=1:size(Team_Names{n}{i},2)%should always be 2
             %[~,~,~,strings]=regexp(Team_Names{n}{i}{j},'[A-Z][a-z][a-z]');
             teamstring=Team_Names{n}{i}{1};
             T_Codes{n,1}=teamstring(1:3);
             T_Codes{n,2}=teamstring(4:end);
             
         %end
     end
 end




counter=1;
for boxscore=1:size(Stats,2)
    
    for i=1:4
        
        for row=2:size(Stats{boxscore}{i},1)
            
            if  isequal((Stats{1,boxscore}{1,i}(row,1)),{'Starters'})||isequal((Stats{1,boxscore}{1,i}(row,1)),{'Reserves'})||isequal((Stats{1,boxscore}{1,i}(row,1)),{'Team Totals'})
                   Stats{1,boxscore}{1,i}(row,1)
                   display('skipping line');
                   

               
            else
                   if i<3 full{counter,1}=T_Codes{boxscore,1}; full{counter,2}=T_Codes{boxscore,2}; 
                       
                   else full{counter,1}=T_Codes{boxscore,2}; full{counter,2}=T_Codes{boxscore,1};end
                 
                   for col=1:size(Stats{boxscore}{i},2)
               
                
                    
                        if ~isequal((Stats{1,boxscore}{1,i}(row,col)),{[]});
                        
                        full{counter,col+2}=cell2mat(Stats{1,boxscore}{1,i}(row,col));
                    
                        else
                        
                        full{counter,col+2}=cell2mat({'NaN'});
                        
                        end
                 end
                 counter=counter+1; 
            end
           
        end
            
        end
    end
matdata=full;
%matdata=cellfun(@cell2mat,full,'UniformOutput',false);
%names={matdata{:,2}};
for i=1:size(matdata,1) names{i,1}=char(matdata{i,3}); end
[C,ic,ia]=unique(names);

%C=cellfun(@cell2mat,C);

new_data=cell(size(ic,1),1);
for i=1:size(ia,1)
    if isempty(new_data{ia(i)})
        
    
        new_data{ia(i)}={char(C{ia(i)}) matdata{i,1} matdata{i,2} char(fake_date{day_n})  matdata{i,4:end}};
    
    else
        new_data{ia(i)}={new_data{ia(i)} matdata{i,5:end}};
    end
end
j=1;
for n=1:size(new_data,1) 
    for i=1:size(new_data{n}{1},2)

                    
                finaldata{n,i}=new_data{n}{1}{i};
                
              
           end
end
    for n=1:size(new_data,1)
        for i=2:size(new_data{n},2)
            finaldata{n,i+size(new_data{n}{1},2)}=new_data{n}{i};
           
        end
    end

savestring=sprintf('%s_fbb.csv',fake_date{day_n});
cell2csv(savestring,finaldata);
end
end
% for k=1:size(Teams)
%     for i=1:size(Teams{k})
%         for j=1:size(daily_stats,2)
%          [x y c]=find(strcmp(Teams{k}{i},Stats{1,1}{1,1})==1);
%          
%          if ~isempty(x)
%              player_loc{i}=[x,1,j];
%          else
%          [x y c]=find(strcmp(Teams{k}{i},Stats{1,1}{1,3})==1);
%          if ~isempty(x)
%              player_loc{i}=[x,3,j];
%          end
%          
%          
%          
%          end
%          
%         end
%         if ~isempty(player_loc)
%         
%         
%         end
%     end
% end