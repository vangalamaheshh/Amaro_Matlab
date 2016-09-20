function draw_total_mutation_rate_graph(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'include_silent_in_sort',true);
P=impose_default_value(P,'exclude_these_genes_from_silent_calculation',{});
P=impose_default_value(P,'exclude_top_n_genes',40);
P=impose_default_value(P,'manual_correction_TSP_WashU',false);
P=impose_default_value(P,'threshold_line',11);
P=impose_default_value(P,'filename','total_rate_graph.jpg');
P=impose_default_value(P,'dimensions',[833 602]);
P=impose_default_value(P,'ymax',0.7e-5);
P=impose_default_value(P,'legend_x',10);

C = compute_mutation_rates(M,P);

figure(1);
close(1);
figure(1);
first_g = find(C.Rnon(C.nc,C.HAT,:)>0|C.Rsil(C.nc,C.HAT,:)>0,1);
ymax = P.ymax;

ylabel('mutation rate');
xlabel('number of top genes omitted');
set(gca,'XLim',[0 M.ng-first_g]);
set(gca,'XDir','reverse');
set(gca,'YLim',[0 ymax]);

legx=P.legend_x;
text(M.ng-first_g-legx,ymax*0.95,'Silent','color',[0 0 1]);
text(M.ng-first_g-legx,ymax*0.90,'Nonsilent','color',[0 0 0]);

hold on
for g=1:M.ng
  c=C.nc;  % total
  x=M.ng-g;

  %nonsilent
  line([x-0.5 x+0.5],[C.Rnon(c,C.HAT,g) C.Rnon(c,C.HAT,g)],'color',[0 0 0]);
  if C.Rnon(c,C.HAT,g)>0
    line([x x],[C.Rnon(c,C.LOW,g) C.Rnon(c,C.HIGH,g)],'color',[0 0 0]);
  end

  %silent
  line([x-0.5 x+0.5],[C.Rsil(c,C.HAT,g) C.Rsil(c,C.HAT,g)],'color',[0 0 1]);
  if C.Rsil(c,C.HAT,g)>0
    line([x x],[C.Rsil(c,C.LOW,g) C.Rsil(c,C.HIGH,g)],'color',[0 0 1]);
  end

end

threshold = P.threshold_line;

line([1 1]*(threshold-0.5),[0.02 0.9]*ymax,'color',[1 0 0]);

p = get(gcf,'position');
pp = get(gcf, 'paperposition');
res = p(3)/pp(3);
w = P.dimensions(1);
h = P.dimensions(2);

set(gcf,'position',[405 78 w h],'color',[1 1 1]);
set(gcf,'paperposition',[0 0 w h]/res);
set(gca,'position',[0.05 0.11 0.9 0.815]);
set(gcf,'papersize',[w h]/res);
hold off

%% PRINT TO FILE

[dev res] = interpret_print_filename(P.filename);
fprintf('Outputting histogram to %s\n', P.filename);
print(['-d' num2str(dev)], ['-r' num2str(res)], P.filename);
