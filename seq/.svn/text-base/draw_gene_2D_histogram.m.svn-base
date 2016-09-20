function draw_gene_2D_histogram(M,P)

% default parameters

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'hist_with_error_bars',false);
P=impose_default_value(P,'hist_with_trend_lines',false);
P=impose_default_value(P,'genes_to_color',{});
P=impose_default_value(P,'gene_histogram_filename','gene_histogram.png');
P=impose_default_value(P,'axismax',500);
P=impose_default_value(P,'label_position_manual_adjustment',{});

fprintf('Drawing mutation-rate 2D histogram\n');

% calculate total coverage for each gene
% and break down into silent and nonsilent "bases at risk"

N_breakdown = repmat(M.N_cov(:,M.TOT,:),1,24) .* repmat(M.breakdown.frac, [1 1 M.np]);
N_si = round(squeeze(sum(sum(N_breakdown(:,1:8,:),2),3)));
N_ns = round(squeeze(sum(sum(N_breakdown(:,9:24,:),2),3)));

n_ns = squeeze(sum(M.n_nonsilent(:,M.TOT,:),3));
n_si = squeeze(sum(M.n_silent(:,M.TOT,:),3));

% compute rates and confidence intervals

[phat_ns cdf_ns] = binofit(n_ns, N_ns, 0.05);
[phat_si cdf_si] = binofit(n_si, N_si, 0.05);

% generate figure

close all;
hold on

axismax = P.axismax;

sc=1000000;

% make sure data falls within limits
if max(phat_si)*sc > axismax || max(phat_ns)*sc > axismax
   fprintf('Warning: Some histogram data falls outside axismax=%d\n', axismax);
end

% set(gca,'YDir','reverse');
ylabel('nonsilent mutations / MB-at-risk');
xlabel('silent mutations / MB-at-risk');
set(gca,'XLim',[0 axismax]);
set(gca,'YLim',[0 axismax]);

% find fits

x=phat_si;y=phat_ns;a=1;
b = [-40:0.01:40]/sc;
ssd = zeros(length(b),1);
for i=1:length(b);
  y_model = (a*x)+b(i);
  sd = (y - y_model) .^ 2;
  ssd(i) = sum(sd);
end
b_best = b(find(ssd==min(ssd)));

x=cdf_si(:,1);y=cdf_ns(:,1);a=1;
b = [-40:0.01:40]/sc;
ssd = zeros(length(b),1);
for i=1:length(b);
  y_model = (a*x)+b(i);
  sd = (y - y_model) .^ 2;
  ssd(i) = sum(sd);
end
b_low = b(find(ssd==min(ssd)));

x=cdf_si(:,2);y=cdf_ns(:,2);a=1;
b = [-40:0.01:40]/sc;
ssd = zeros(length(b),1);
for i=1:length(b);
  y_model = (a*x)+b(i);
  sd = (y - y_model) .^ 2;
  ssd(i) = sum(sd);
end
b_high = b(find(ssd==min(ssd)));

if P.hist_with_trend_lines
  line([0 axismax]*sc, ([0 axismax]+b_low)*sc,'Color',[0 0 0]);
  line([0 axismax]*sc, ([0 axismax]+b_best)*sc,'Color',[0 1 0]);
  line([0 axismax]*sc, ([0 axismax]+b_high)*sc,'Color',[0 0 0]);
else
  % just draw y=x line
  line([0 axismax]*sc, ([0 axismax]+0)*sc,'Color',[0 1 0]);
end

% plot data

for i=1:M.ng
  color = [0 0 0];
  if ismember(M.gene.name{i},P.genes_to_color), color = [1 0 0]; end

  xcen = phat_si(i);
  xlow = cdf_si(i,1);
  xhi  = cdf_si(i,2);
  ycen = phat_ns(i);
  ylow = cdf_ns(i,1);
  yhi  = cdf_ns(i,2);

  if P.hist_with_error_bars
    line([xlow xhi]*sc,[ycen ycen]*sc,'Color',color);
    line([xcen xcen]*sc,[ylow yhi]*sc,'Color',color);
  else
    r=axismax/120;
    x=(xcen*sc)-r;
    y=(ycen*sc)-r;
    x=round(x*10)/10;
    y=round(y*10)/10;
    rectangle('position',[x y 2*r 2*r],'curvature',[1 1],...
      'EdgeColor', color, 'FaceColor', color, 'clipping', 'off');
  end

end % next point

set_all_fontsize(gcf, 16);
p = get(gcf,'position');
pp = get(gcf, 'paperposition');
res = p(3)/pp(3);
w = 782;
h = 672;
set(gcf,'position',[77 297 w h],'color',[1 1 1]);
set(gcf,'paperposition',[0 0 w h]/res);
set(gcf,'papersize',[w h]/res);
hold off

filename = P.gene_histogram_filename;
[dev res] = interpret_print_filename(filename);
fprintf('Outputting histogram to %s\n', filename);
print(['-d' dev], ['-r' num2str(res)], filename);

end



