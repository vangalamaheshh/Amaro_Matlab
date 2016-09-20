clear p_pick_winner
clear p_pick_loser
v = .1;
n=1;
w=.9;
i=1;
p_pick_loser(1) = v/(v+((n-i+1)*w));
p_pick_winner(1) = ((n-i+1)*w)/(v+((n-i+1)*w));
for i=2:n+1
    p_pick_winner(i) = ((n-i+1)*w)/(v+((n-i+1)*w));
    p_pick_loser(i) = v/(v+(n-i+1)*w) * prod(p_pick_winner(1:i-1));
end



p=[0:.01:1];
d=[0:.01:1];
for i =1:length(d)
y(i,:)=d(i)*p+(1-d(i))*(1-p);
end

