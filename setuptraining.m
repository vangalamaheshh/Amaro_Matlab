usertrain.actual_play=zeros(slength(usertrain),1);
order_m=unique(profiles.user);

for i=1:slength(usertrain)
    l=find(ismember(order_m,usertrain.user{i}));
    vals=find(m(l,:)>0);
    x=randi(length(vals),1);
    usertrain.actual_play(i,1)=m(l,vals(x));
    if mod(i,10000)==0
        i
    end
end