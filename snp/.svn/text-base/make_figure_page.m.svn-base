function make_figure_page(perpage_x,perpage_y,Ch,Cc,i,maxii)
close all
f1=figure(1);
set(gcf,'Visible','off');
gr=make_subplotgrid(4*ones(1,perpage_x),3*ones(1,perpage_y),...
                    0.5*ones(1,perpage_x+1),1.2*ones(1,perpage_y+1),0.3,0.45);

for ii=i:min(i+perpage_y-1,maxii)
  disp(ii);
  
  a=subplotgrid(gr,ii-i+1,1);
  h=hgload(['Sample' sprintf('%03d',ii) '.smoothed' '.fig']);
  set(gcf,'Visible','off');    
  ylabel('log2ratio','FontSize',3);
  fo=findobj(h,'Type','Axes');
  c=copyobj(fo,f1);
  set(c,'Position',get(a,'Position'));
  delete(a);
  axes(c);
  set(gcf,'Visible','off');
  colormap(1-0.5*repmat((0:(1/63):1)',1,3));
  close(h);
  
  a=subplotgrid(gr,ii-i+1,2);
  h=hgload(['Sample' sprintf('%03d',ii) '.hist' '.fig']);
  set(gcf,'Visible','off');    
  xlabel('log2ratio','FontSize',3);
  smc=sort([ Ch.peaks{ii} Cc.peaks{ii}  ]);
  ax=axis;
  ax=[ smc(1)-(smc(end)-smc(1))*0.1-0.1 smc(end)+(smc(end)-smc(1))*0.1+0.1 0 ax(4)];
  axis(ax);
  if (0)
    fo=findobj(h,'Type','line');
    fo1=findobj(fo,'Color',[1 0 0]);
    for kk=1:length(fo1)
      yd=get(fo1(kk),'YData');
      set(fo1(kk),'YData',[0 yd(2)*0.2]);
      ax(4)=yd(2);
    end
    fo1=findobj(fo,'Color',[0 1 0]);
    for kk=1:length(fo1)
        yd=get(fo1(kk),'YData');
        set(fo1(kk),'YData',[yd(2)*0.25 yd(2)*0.45]);
    end    
    axis(ax);
  end
  fo=findobj(h,'Type','Axes');    
  c=copyobj(fo,f1);
  set(c,'Position',get(a,'Position'));
  delete(a);
  axes(c);
  set(gcf,'Visible','off');
  close(h);
  
  a=subplotgrid(gr,ii-i+1,3);
  h=hgload(['Sample' sprintf('%03d',ii) '.levels' '.fig']);
  set(gcf,'Visible','off');    
  ylabel('log2ratio','FontSize',3);
  set(gca,'LineWidth',0.5);
  set(gca,'FontSize',5);
  fo=findobj(gcf,'Type','Line','Color',[1 0 0]);
  set(fo,'LineWidth',0.5);
  fo=findobj(gcf,'Type','Line','Color',[0 1 0]);
  set(fo,'LineWidth',0.4);
  fo=findobj(h,'Type','Axes');    
  c=copyobj(fo,f1);
  set(c,'Position',get(a,'Position'));
  delete(a);
  axes(c);
  set(gcf,'Visible','off');        
  colormap(1-0.5*repmat((0:(1/63):1)',1,3));
  close(h);
  
  if (0)
    a=subplotgrid(gr,ii-i+1,4);
    h=hgload(['Sample' sprintf('%03d',ii) '.loh' '.fig']);
    c=copyobj(get(h,'Children'),f1);
    set(c,'Position',get(a,'Position'));
    delete(a);
    axes(c);
    close(h); 
  end
  
  a=subplotgrid(gr,ii-i+1,4);
  h=hgload(['Sample' sprintf('%03d',ii) '.barpie' '.fig']);
  set(gcf,'Visible','off');    
  xlabel('log2ratio','FontSize',3);
  set(gca,'LineWidth',0.1);
  set(gca,'FontSize',5);
  cax=axis;
  axis([ax(1:2) cax(3:4)]);
  fo=findobj(gcf,'Type','Text');
  set(fo,'FontSize',4);
  fo=findobj(gcf,'Type','Line','YData',[ 0 0]);
  set(fo,'Xdata',ax(1:2));
  fo=findobj(h,'Type','Axes');    
  c=copyobj(fo,f1);
  set(c,'Position',get(a,'Position'));
  delete(a);
  axes(c);
  set(gcf,'Visible','off');        
  close(h); 
  
  a=subplotgrid(gr,ii-i+1,5);
  h=hgload(['Sample' sprintf('%03d',ii) '.cn' '.fig']);
  xlabel('2*ratio','FontSize',3);
  ylabel('corrected copy #','FontSize',3);
  lh=legend;
  set(gca,'LineWidth',0.1);
  set(gca,'FontSize',5);
  delete(lh);
  fo=findobj(gcf,'Type','line','Color',[0.8 0.8 0.8]);
  cax=axis;
  %      axis([ get(fo(1),'XData') floor(cax(3))-0.1 ceil(cax(4))+0.1]);
  xd=get(fo(1),'XData');
  axis([ xd(1) min(xd(2),5) floor(cax(3))-0.1 min(ceil(cax(4))+0.1,5)]);
  fo=findobj(gcf,'Type','line');
  set(fo,'MarkerSize',3);
  fo=findobj(gcf,'Type','Text');
  set(fo,'FontSize',4);
  fo=findobj(h,'Type','Axes');    
  c=copyobj(fo,f1);
  set(c,'Position',get(a,'Position'));
  delete(a);
  axes(c);
  set(gcf,'Visible','off');
  close(h); 
end
