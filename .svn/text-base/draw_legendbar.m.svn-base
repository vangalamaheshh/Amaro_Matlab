function draw_legendbar(horiz,convtab,posvec,textvec,fs)

cmat=zeros(size(convtab,1),1,3);
for i=1:size(convtab,1)
  cmat(i,1,:)=convtab{i,2};
end

if strcmp(horiz(1:4),'hori')
    pos=get(gca,'Position');
    pos2=[pos(1) pos(2)+pos(4)/2 pos(3) pos(4)/2];
    set(gca,'Position',pos2);
    image(permute(cmat,[2 1 3])); 
    for i=1:length(posvec)
        text(posvec(i),2,textvec{i},'FontSize',fs,'HorizontalAlignment','center');
    end 
    set(gca,'XTick',posvec,'XTickLabel',[]);
    set(gca,'YTick',[],'XTickLabel',[]);    
    set(gca,'TickLength',[0 0]);
else
    image(cmat);
    set(gca,'YDir','normal');
    set(gca,'YTick',posvec,'YTickLabel',textvec,'FontSize',fs,'XTick',[]);
    set(gca,'TickLength',[0 0]);
end
