function [interval,pos]=map_to_intervals(idx,start,bpi)

interval = zeros(length(idx),1);
pos = zeros(length(idx),1);
bpi = [0; bpi];
for i=1:length(idx)
  interval(i) = find(bpi>=idx(i),1)-1;
  pos(i) = idx(i)+start(interval(i))-bpi(interval(i))-1;
end

if 0
  both=zeros(length(idx)+length(bpi),2);
  both(:,1)=[as_column(idx); as_column(bpi)];
  both(length(idx)+1:end,2)=1;
  
  [tmp,si]=sort(both(:,1));
  %% cumsum ... 
end
