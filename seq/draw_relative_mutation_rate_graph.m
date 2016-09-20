function draw_relative_mutation_rate_graph(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'include_silent_in_sort',true);
P=impose_default_value(P,'exclude_these_genes_from_silent_calculation',{});
P=impose_default_value(P,'exclude_top_n_genes',40);
P=impose_default_value(P,'manual_correction_TSP_WashU',false);
P=impose_default_value(P,'threshold_line',11);
P=impose_default_value(P,'filename','total_rate_graph.png');
P=impose_default_value(P,'dimensions',[824 618]);
P=impose_default_value(P,'legend_x',20);
P=impose_default_value(P,'ymax',12);

C = compute_mutation_rates(M,P);

% set up graph
figure(2);
close(2);
figure(2);
first_g = find(~isnan(C.RRnon(1,C.HAT,:)),1);
ymax = P.ymax;
hold on

ylabel('mutation rate (relative to total)');
xlabel('number of top genes omitted');
set(gcf,'color',[1 1 1]);
set(gca,'XLim',[0 M.ng-first_g]);
set(gca,'XDir','reverse');
set(gca,'YLim',[0 ymax]);

black=[0 0 0];

if C.nc-1==3
  colors = { [0 0 1] [0 0.5 0] [1 0 0] };
elseif C.nc-1==5
  colors = { [0 0 1] [0.25 0.6 0] [0 0.6 0.25] [0.8 0.25 0] [0.8 0 0.25] };
elseif C.nc-1==4 && ~isempty(grep(M.mutclass,'C+G',1))
  colors = { [0 0 1] [0 0.5 0] [0.8 0.25 0] [0.8 0 0.25] };
elseif C.nc-1==4 && ~isempty(grep(M.mutclass,'A+T',1))
  colors = { [0 0 1] [0.25 0.6 0] [0 0.6 0.25] [1 0 0] };
else
  for c=1:C.nc-1
    colors{c} = rand(1,3);
  end
end

% add legend
legx = P.legend_x;
for c=1:C.nc-1
  text(M.ng-first_g-legx,ymax-0.3*c,M.mutclass{c},'color',colors{c});
end

% plot nonsilent data

stagger=0.1;
for g=1:M.ng
  for c=1:C.nc-1
    x=M.ng-g+(c*stagger);
    line([x-0.5 x+0.5],[C.RRnon(c,C.HAT,g) C.RRnon(c,C.HAT,g)],'color',colors{c});
    line([x x],[C.RRnon(c,C.LOW,g) C.RRnon(c,C.HIGH,g)],'color',colors{c});
  end
end

% plot silent data
g=M.ng+10;
text(M.ng-g+3,-0.2,'silent','color',black);
for c=1:C.nc-1
  x=M.ng-g+(c*stagger);
  line([x-0.5 x+0.5],[C.RRTsil(c,C.HAT) C.RRTsil(c,C.HAT)],'color',colors{c},'clipping','off');
  line([x x],[C.RRTsil(c,C.LOW) C.RRTsil(c,C.HIGH)],'color',colors{c},'clipping','off');
end

p = get(gcf,'position');
pp = get(gcf, 'paperposition');
res = p(3)/pp(3);
w = P.dimensions(1);
h = P.dimensions(2);

set(gcf,'position',[400 400 w h],'color',[1 1 1]);
set(gcf,'paperposition',[0 0 w h]/res);
set(gcf,'papersize',[w h]/res);
hold off

%% PRINT TO FILE

[dev res] = interpret_print_filename(P.filename);
fprintf('Outputting histogram to %s\n', P.filename);
print(['-d' num2str(dev)], ['-r' num2str(res)], P.filename);

