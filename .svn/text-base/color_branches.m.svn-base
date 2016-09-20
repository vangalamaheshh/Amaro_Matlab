function color_branches(lhs,lnk,labels,c)

u=unique(labels);
for i=1:length(u)
    mn=min(find(labels==u(i)));
    mx=max(find(labels==u(i)));
    ymax=max(lnk(find(all(lnk(:,1:4)<=mx & lnk(:,1:4)>=mn,2)),5));
%    disp([ u(i) mn mx ymax]);
    for j=1:length(lhs)
        xd=get(lhs(j),'XData');
        yd=get(lhs(j),'YData');
        if all(xd<=mx & xd>=mn) && all(yd<=ymax)
%            disp([ xd yd ]);
            set(lhs(j),'Color',c(u(i),:));
        end
    end
%    pause
end 
