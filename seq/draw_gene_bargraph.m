function draw_gene_bargraph(M,P)

% default parameters

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'bargraph_genes',{});
P=impose_default_value(P,'bargraph_filename','gene_bargraph.png');
P=impose_default_value(P,'bargraph_fontsize',16);
P=impose_default_value(P,'bargraph_grey_horiz_rules',true);
P=impose_default_value(P,'bargraph_display_mutation_table',false);
P=impose_default_value(P,'bargraph_display_counts_in_bars',false);
P=impose_default_value(P,'bargraph_display_counts_to_right',false);
P=impose_default_value(P,'bargraph_white_numbers',false);

fprintf('Drawing gene bargraph\n');

gene = M.mut.gene(M.use_nonsilent);
patient = M.mut.patient(M.use_nonsilent);
is_treated = M.patient.treated(patient);
is_hypermutated = M.patient.hypermutated(patient);

genes = listmap(P.bargraph_genes,M.gene.name);
genes = genes(~isnan(genes));
ngenes = length(genes);

NON=1;  % non-treated, non-hypermutated
TRE=2;  % treated, non-hypermutated
HYP=3;  % hypermutated

n = zeros(ngenes,3);

for g=1:ngenes
  n(g,NON) = sum(gene==genes(g) & ~is_treated & ~is_hypermutated);
  n(g,TRE) = sum(gene==genes(g) & is_treated & ~is_hypermutated);
  n(g,HYP) = sum(gene==genes(g) & is_hypermutated);
end

nmax = max(sum(n,2));
ymax = 5*ceil((nmax+2)/5);

close all;
hold on

set(gca,'xlim',[0 ngenes]);
set(gca,'ylim',[0 ymax]);
set(gca,'xticklabel',{});
set(gca,'tickdir','out');

if P.bargraph_display_mutation_table
  set(gca,'position',[0.2 0.22 0.75 0.72]);
  h = 600;
  w = 1000;
  texty0 = -0.7*ymax/15;
  textysp = 0.8*ymax/15;
  labelx = -1;
else
  set(gca,'position',[0.13 0.15 0.8 0.8]);
  h = 480;
  w = 850;
  texty0 = -0.8*ymax/15;
  textysp = 1*ymax/15;
  labelx = -0.5;
end

p = get(gcf,'position');
pp = get(gcf, 'paperposition');
res = p(3)/pp(3);
set(gcf,'position',[0 0 w h],'color',[1 1 1]);
set(gcf,'paperposition',[0 0 w h]/res);
set(gcf,'papersize',[w h]/res);

set_all_fontsize(gcf, P.bargraph_fontsize);

qy = texty0 - textysp;
text(labelx,qy,'q value:','HorizontalAlignment','Center',...
     'fontsize',P.bargraph_fontsize-2);

purple = [0.45 0 0.45; 0.7 0.3 0.7; 0.83 0.63 0.93];
categ = {'Non-treated','Treated','Treated + Hypermutated'};

if P.bargraph_display_mutation_table
    for i=1:3
      y = qy - i*textysp;
      text(labelx,y,categ{4-i},'HorizontalAlignment','Center',...
       'fontsize',P.bargraph_fontsize-2);
    end
end  

if P.bargraph_grey_horiz_rules
  for y=5:5:ymax
    line([0 ngenes],[y y],'color',[0.8 0.8 0.8]);
  end
end

set(gca,'ytick',[0:5:ymax]);

text(-0.65,ymax/2,'Somatic Mutations','rotation',90, 'fontsize', P.bargraph_fontsize, ...
   'horizontalalignment','center','verticalalignment','middle');

barwidth = 0.38;

for g=1:ngenes

  % draw stacked bar

  y1 = n(g,NON);
  y2 = y1+n(g,TRE);
  y3 = y2+n(g,HYP);
  x1 = g-0.5-0.5*barwidth;

  if y1, rectangle('position',[x1 0 barwidth y1], 'facecolor', purple(1,:)); end
  if y2-y1, rectangle('position',[x1 y1 barwidth (y2-y1)], 'facecolor', purple(2,:)); end
  if y3-y2, rectangle('position',[x1 y2 barwidth (y3-y2)], 'facecolor', purple(3,:)); end

  % write gene name

  text(g-0.5,texty0,P.bargraph_genes{g},'HorizontalAlignment','Center',...
      'fontangle','italic','fontsize',P.bargraph_fontsize);

  % write q value

  q = M.Q(genes(g));
  if q==0
    Q_text = '<10^-^8';
  else
    exp = -floor(log10(q));
    man = q*10^exp;
    Q_text = sprintf('%1.0fx10^-^%d',man,exp);
  end
  text(g-0.5,qy,Q_text,'HorizontalAlignment','Center',...
     'fontsize',P.bargraph_fontsize-2);

  if P.bargraph_display_mutation_table
    for i=1:3
      y = qy - i*textysp;
      text(g-0.5,y,num2str(n(g,4-i)),'HorizontalAlignment','Center',...
       'fontsize',P.bargraph_fontsize-2);
    end
  end

  if P.bargraph_display_counts_in_bars
    for i=1:3
     ct = n(g,4-i);
     if ~ct, continue; end
     color = [0 0 0];
     if i==1
       if ct<2, y = y3 + 1;
       else y = y2+(y3-y2)/2; end
     elseif i==2
       if ct<2, y = y2 + 1;
       else y = y1+(y2-y1)/2; end
     else
       y = y1/2;
       if P.bargraph_white_numbers, color = [1 1 1]; end
     end
     text(g-0.5,y,num2str(ct),'HorizontalAlignment','Center',...
       'fontsize',P.bargraph_fontsize-2,'color',color);
    end
  end

  if P.bargraph_display_counts_to_right
    for i=1:3
     ct = n(g,4-i);
     if ~ct, continue; end
     color = [0 0 0];
     if i==1
       y = y2+(y3-y2)/2;
     elseif i==2
       y = y1+(y2-y1)/2;
     else
       y = y1/2;
     end
     text(g-0.25,y,num2str(ct),'HorizontalAlignment','left',...
       'fontsize',P.bargraph_fontsize-6,'color',color);
    end
  end

end

% draw legend

lx = 4.5 * (ngenes/8);
ly = 28 * (ymax/40);
lw = 3.1 * (ngenes/8);
lh = 10 * (ymax/40);
sw = 0.2 * (ngenes/8);
sh = 1.9 * (ymax/40);

rectangle('position',[lx ly lw lh],'facecolor',[1 1 1],'edgecolor',[0 0 0]);

for y=1:3
  yy = ly+(lh*0.05)+(y-0.85)*((lh*0.9)/3); 
  rectangle('position', [lx+0.2, yy, sw, sh],...
     'edgecolor',[0 0 0],'facecolor',purple(y,:));
  text(lx+0.4+sw, yy, categ{y},'verticalalignment','bottom',...
     'fontsize',P.bargraph_fontsize);
end

% output figure

hold off

[dev res] = interpret_print_filename(P.bargraph_filename);
fprintf('Outputting histogram to %s\n', P.bargraph_filename);
print(['-d' dev], ['-r' num2str(res)], P.bargraph_filename);

end
