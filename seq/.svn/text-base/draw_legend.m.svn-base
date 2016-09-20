function success = draw_legend(P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'cluster_labels','OV');   % GBM or OV

fontsize = 5;
% SET UP PLOT AREA

figure(1)
clf
set(gcf,'visible','off');
%set(gcf,'visible','on');
set(gcf,'renderer','zbuffer');
set(gcf,'renderermode','auto');
plot(0,0);    % needed for top and right black lines
axis([0 1 0 1]);
set(gca,'position',[0 0 1 1]);
set(gca,'xticklabel',{});
set(gca,'yticklabel',{});
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontSize',fontsize);
hold on

% EDGE COLOR = mutation

white = [1 1 1];
black = [0 0 0];
dkgrey = 0.5*white;
ltgrey = 0.8*white;
red = [1 0 0];
green = [0.3 0.8 0];
blue = [0 0 1];
orange = [0.9 0.5 0];
yellow = [0.9 0.8 0];
pink = [1 0.5 0.6];
purple = [0.5 0 0.5];

truncated_mut_color = red;
other_mut_color = pink;
methylated_scheme = [0.9 0.7 0;0 0.3 0;0 0 0.8];

% POSITIONS OF LEGEND ITEMS

x1 = 0.03;
x2 = 0.10;
y = 0.97;
yincr = -0.038;

% DRAW LEGEND

% methylation colors

draw_shape('cir',x1,y,[1 1 0]*methylated_scheme);
text(x2,y,'Fully methylated');
y=y+yincr;

draw_shape('cir',x1,y,[0.5 1 0.5]*methylated_scheme);
text(x2,y,'Intermed. methylation');
y=y+yincr;

draw_shape('cir',x1,y,[0 1 1]*methylated_scheme);
text(x2,y,'Unmethylated');
y=y+yincr;

draw_shape('cir',x1,y,dkgrey);
text(x2,y,'No methylation data');
y=y+yincr;

y=y+0.5*yincr;

% mutation halos

draw_shape('cir',x1,y,white,'small',truncated_mut_color);
text(x2,y,'Truncating mutation');
y=y+yincr;

draw_shape('cir',x1,y,white,'small',other_mut_color);
text(x2,y,'Other mutation');
y=y+1.15*yincr;

draw_shape('cir',x1,y,white,'large',truncated_mut_color);
text(x2,y,'2+ muts (1+ trunc.)');
y=y+1.3*yincr;

draw_shape('cir',x1,y,white,'large',other_mut_color);
text(x2,y,'2+ muts (none trunc.)');
y=y+1.15*yincr;

draw_shape('cir',x1,y,white,'small',ltgrey);
text(x2,y,'No mutations');
y=y+yincr;

rectangle('position',[x1-0.014/2 y-0.014/2 0.014 0.014],'facecolor',white,'curvature',[1 1]);
rectangle('position',[x1-0.026/2 y-0.026/2 0.026 0.026],'facecolor',white,'curvature',[1 1]);
text(x2,y,'Not sequenced');
y=y+yincr;

y=y+0.5*yincr;

% treated/hypermutated status

draw_shape('cir',x1,y,dkgrey);
draw_shape('dot',x1,y,black);
text(x2,y,'Treated sample');
y=y+yincr;

draw_shape('cir',x1,y,dkgrey);
draw_shape('dot',x1,y,red);
text(x2,y,'Hypermutated sample');
y=y+yincr;


y=y+0.5*yincr;

% cluster shapes

switch(P.cluster_labels)
  case 'GBM', clabs = {'Proneural cluster','Neural cluster','Classical cluster',...
                        'Mesenchymal cluster','Cluster unknown'};
  case 'OV',  clabs = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster unknown'};
  otherwise, error('Unknown cluster');
end

draw_shape('sqr',x1,y,black);
text(x2,y,clabs{1});
y=y+yincr;

draw_shape('inv',x1,y,black);
text(x2,y,clabs{2});
y=y+yincr;

draw_shape('tri',x1,y,black);
text(x2,y,clabs{3});
y=y+yincr;

draw_shape('dmd',x1,y,black);
text(x2,y,clabs{4});
y=y+yincr;

draw_shape('cir',x1,y,black);
text(x2,y,clabs{5});
y=y+yincr;




% OUTPUT LEGEND

set_all_fontsize(gca,fontsize);
set(gcf,'paperposition',[0 0 2.4 2.5]);
hold off
print('-dpng','-r250','-zbuffer','legend.png');
success = true;
return

