function venn(s123,s100,s200,s300,s120,s023,s103,s000)
% venn(s123,s100,s200,s300,s120,s023,s103,s000)
%
% based on venn_diagram.m (Gaddy Getz)

  colors = {[0 0 1],[0.2 0.7 0.2],[1 0 0]};
  places = {[0.12 0.25],[0 0],[0.25 0]};

  clf,hold on

  r1=rectangle('Position',[places{1} 0.5 0.5],'Curvature',[1 1],'linewidth',2, ...
               'EdgeColor',colors{1});
  r2=rectangle('Position',[places{2} 0.5 0.5],'Curvature',[1 1],'linewidth',2, ...
               'EdgeColor',colors{2});
  r3=rectangle('Position',[places{3} 0.5 0.5],'Curvature',[1 1],'linewidth',2, ...
               'EdgeColor',colors{3});

  pos=get(r1,'Position'); 
  plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
       'Color',get(r1,'EdgeColor'));
  pos=get(r2,'Position'); 
  plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
       'Color',get(r2,'EdgeColor'));
  pos=get(r3,'Position'); 
  plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
       'Color',get(r3,'EdgeColor'));
  
  axis square
  axis off

  arg = {'fontsize',14,'horizontalalignment','center'};

  text(0.35,0.65,sub_num2str(s100),arg{:});
  text(0.20,0.40,sub_num2str(s120),arg{:});
  text(0.12,0.20,sub_num2str(s200),arg{:});
  text(0.50,0.40,sub_num2str(s103),arg{:});
  text(0.35,0.35,sub_num2str(s123),arg{:});
  text(0.37,0.15,sub_num2str(s023),arg{:});
  text(0.57,0.20,sub_num2str(s300),arg{:});
 
if exist('s000','var')
  tot=s000+s123+s100+s200+s300+s103+s120+s023;
  text(0.00,0.60,{'tot',sub_num2str(tot)},arg{:});
  text(0.75,0.00,{'none',sub_num2str(s000)},arg{:});
end


  hold off

  function txt = sub_num2str(num)
    txt = format_number(num,2,5);
  end



end
