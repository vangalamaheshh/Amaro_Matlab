function [ph,lh]=cut_bar(h,dh,th,gh,col)

b=findobj(gca,'Type','hggroup');
ax=axis;
axis([ax(1:2) ax(3) h+dh]);

for i=1:length(b)
  set(b(i),'Clipping','off');
end
ph=[];
kh=[];

if (0)
  ph=[];
  for i=1:length(b)
    x=get(b(i),'XData');
    y=get(b(i),'YData');
    w=get(b(i),'BarWidth');
    for j=1:length(x)
      if y(j)>(h+th+gh)
        ph(end+1)=patch([x(j)-w/2 x(j)+w/2 x(j)+w/2 x(j)-w/2 x(j)-w/2 ],[h h+th h+th+gh h+gh h],col);
        set(ph(end),'EdgeColor',col,'FaceColor',col,'LineStyle','none');
        if (0)
          line([x(j)-w/2 x(j)+w/2],[h h+th],'LineWidth',get(b(i),'LineWidth'),'LineStyle',get(b(i),'LineStyle'),...
               'Color',get(b(i),'EdgeColor'));
          line([x(j)-w/2 x(j)+w/2],[h+gh h+th+gh],'LineWidth',get(b(i),'LineWidth'),'LineStyle',get(b(i),'LineStyle'),...
               'Color',get(b(i),'EdgeColor'));
        end
      end
    end
  end
end
