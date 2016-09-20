function sep=venn_diagram(groups,desc)

if length(groups)==2
%    figure;
    r1=rectangle('Position',[0 0 0.5 0.5],'Curvature',[1 1],'linewidth',2,...
        'EdgeColor',[1 0 0]);
    r2=rectangle('Position',[0.25 0 0.5 0.5],'Curvature',[1 1],'linewidth',2,...
        'EdgeColor',[0 0 1]);
    hold on
    pos=get(r1,'Position'); 
    plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
        'Color',get(r1,'EdgeColor'));
    pos=get(r2,'Position'); 
    plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
        'Color',get(r2,'EdgeColor'));
    legend(desc);
    axis square
    axis off
    sep{1}=setdiff(groups{1},groups{2});
    sep{2}=setdiff(groups{2},groups{1});
    sep{3}=intersect(groups{1},groups{2});
    text(0.12,0.25,num2str(length(sep{1})),'FontSize',14);
    text(0.37,0.25,num2str(length(sep{3})),'FontSize',14);
    text(0.67,0.25,num2str(length(sep{2})),'FontSize',14);
elseif length(groups)==3
%    figure;
    r1=rectangle('Position',[0 0 0.5 0.5],'Curvature',[1 1],'linewidth',2, ...
        'EdgeColor',[1 0 0]);
    r2=rectangle('Position',[0.25 0 0.5 0.5],'Curvature',[1 1],'linewidth',2, ...
        'EdgeColor',[0 1 0]);
    r3=rectangle('Position',[0.12 0.25 0.5 0.5],'Curvature',[1 1],'linewidth',2, ...
       'EdgeColor',[0 0 1]);
    hold on
    pos=get(r1,'Position'); 
    plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
        'Color',get(r1,'EdgeColor'));
    pos=get(r2,'Position'); 
    plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
        'Color',get(r2,'EdgeColor'));
    pos=get(r3,'Position'); 
    plot([ 0.5*(pos(1)+pos(3)) 0.5*(pos(1)+pos(3))],[pos(2) pos(2)],...
        'Color',get(r3,'EdgeColor'));
    
    legend(desc);
    axis square
    axis off
    sep{1}=setdiff(groups{1},union(groups{2},groups{3}));
    sep{2}=setdiff(groups{2},union(groups{1},groups{3}));
    sep{3}=setdiff(groups{3},union(groups{1},groups{2}));
    sep{7}=intersect(intersect(groups{1},groups{2}),groups{3});
    sep{4}=setdiff(intersect(groups{1},groups{2}),sep{7});
    sep{5}=setdiff(intersect(groups{2},groups{3}),sep{7});
    sep{6}=setdiff(intersect(groups{3},groups{1}),sep{7});
    
    a1=length(groups{1}); 
    a2=length(groups{2});
    a3=length(groups{3});
    a12=length(intersect(groups{1},groups{2}));
    a13=length(intersect(groups{1},groups{3}));
    a23=length(intersect(groups{2},groups{3}));
    a123=length(intersect(intersect(groups{1},groups{2}),groups{3}));
    
%     a1-a12-a13+a123-length(sep{1})
%     a2-a12-a23+a123-length(sep{2})
%     a12-a123-length(sep{4})
%     a13-a123-length(sep{6})
%     a123-length(sep{7})
%     a23-a123-length(sep{5})
%     a3-a13-a23+a123-length(sep{3})
    
    
    text(0.12,0.2,num2str(a1-a12-a13+a123),'FontSize',14);
    text(0.37,0.15,num2str(a12-a123),'FontSize',14);
    text(0.57,0.2,num2str(a2-a12-a23+a123),'FontSize',14);
    text(0.2,0.4,num2str(a13-a123),'FontSize',14);
    text(0.35,0.35,num2str(a123),'FontSize',14);
    text(0.5,0.4,num2str(a23-a123),'FontSize',14);
    text(0.35,0.65,num2str(a3-a13-a23+a123),'FontSize',14);
    
else
    error('too few or many groups');
end