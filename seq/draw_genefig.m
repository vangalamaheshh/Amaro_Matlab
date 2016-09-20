function [success,maptext] = draw_genefig(M,gene,P)
%
% draws graphic for genepages
%

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'visible','off');
P=impose_default_value(P,'show_unknowns',true);
P=impose_default_value(P,'randseed',1234);
P=impose_default_value(P,'suppress_title',false);
P=impose_default_value(P,'label_top_samples',0);
P=impose_default_value(P,'resolution',250);
P=impose_default_value(P,'expression_axis_label','Normalized expression');

rand('twister',P.randseed);

success = false;
maptext = [];

gidx = find(strcmpi(gene,M.gene.name),1);
if isempty(gidx)
  fprintf('Gene not found: %s\n', gene);
  return
end

fontsize = 5;

% GATHER DATA

D = M.patient;
D.expr = M.expr(gidx,:)';
D.cn = M.cn(gidx,:)';
D.methyl = M.methyl(gidx,:)';
D.have_expr = ~isnan(D.expr);
D.have_cn = ~isnan(D.cn) & ~isinf(D.cn);
D.have_methyl = ~isnan(D.methyl);
D.have_coverage = squeeze(M.N_cov(gidx,M.TOT,:)>0);
D.n_nons = squeeze(M.n_nonsense(gidx,M.TOT,:));
D.n_fshift = squeeze(sum(M.n_indel(gidx,grep('Frameshift',M.mutclass,1),:),2));
D.n_infr = squeeze(sum(M.n_indel(gidx,grep('Inframe',M.mutclass,1),:),2));
D.n_splice = squeeze(M.n_splice(gidx,M.TOT,:));
D.n_mis = squeeze(M.n_missense(gidx,M.TOT,:));
D.n_trunc = D.n_nons + D.n_fshift + D.n_splice;
D.n_other = D.n_mis + D.n_infr;
D.n_tot = D.n_trunc + D.n_other;
D.have_seq = D.have_coverage | D.n_tot;

% SET UP PLOT AREA

figure(1)
clf
set(gcf,'renderer','zbuffer');
set(gcf,'renderermode','auto');
set(gcf,'visible',P.visible);
plot(0,0);    % needed for top and right black lines
axis([0 1 0 1]);

if P.suppress_title
  axispos = [0.14 0.1 0.85 0.89];
  paperpos = [0 0 2.6 2.8];
else
  axispos = [0.14 0.1 0.85 0.82];
  paperpos = [0 0 2.6 3.1];
end

paperres = P.resolution;
set(gca,'position',axispos);
set(gcf,'paperposition',paperpos);
hold on

xunk_size = 0.15;
yunk_size = 0.15;
xscatter = 0.5;
yscatter = 0.5;

if P.show_unknowns && ~all(D.have_cn)
  xunks = true;
  xunk_start = 0;
  xunk_end = xunk_size;
  xkn_start = xunk_size;
  xkn_end = 1; 
  line([xkn_start xkn_start],[0 1],'color',[0 0 0]);
else
  xunks = false;
  xunk_start = 0;
  xunk_end = 0;
  xkn_start = 0;
  xkn_end = 1;
end

if P.show_unknowns && ~all(D.have_expr)
  yunks = true;
  yunk_start = 0;
  yunk_end = yunk_size;
  ykn_start = yunk_size;
  ykn_end = 1;
  line([0 1],[ykn_start ykn_start],'color',[0 0 0]);
else
  yunks = false;
  ynk_start = 0;
  yunk_end = 0;
  ykn_start = 0;
  ykn_end = 1;
end
xkn_size = xkn_end - xkn_start;
ykn_size = ykn_end - ykn_start;

% DETERMINE DATA RANGES

cn_min = -1 - any(D.cn==-2);
cn_max = +1 + any(D.cn==+2);
cn_ncats = cn_max - cn_min + 1;

expr_min = -ceil(max(-D.expr))-0.5;
expr_max = ceil(max(D.expr))+0.5;
expr_range = expr_max - expr_min;

% LABEL AXES

set(gca,'xtick',[xkn_start:(xkn_size/cn_ncats):xkn_end]);
set(gca,'xticklabel',{});
set(gca,'ytick',ykn_start+ykn_size*([expr_min+0.5:1:expr_max-0.5]-expr_min)/expr_range);
set(gca,'yticklabel',{});

xunk_label = 'Unknown';
label_y = -0.04;
if xunks
  text(xunk_end/2,label_y,xunk_label,'fontsize',fontsize,'horizontalalignment','center');
end
cat_labels={'HomoDel','HemiDel','Neutral','LowAmp','HighAmp'};
for cat=cn_min:cn_max
  x = xkn_start + (xkn_size / cn_ncats) * (cat - cn_min + 0.5);
  text(x,label_y,cat_labels{cat+3},'fontsize',fontsize,'horizontalalignment','center');
end

yunk_label = 'Unknown';
label_x = -0.02;
if yunks
  text(label_x,yunk_end/2,yunk_label,'fontsize',fontsize,'verticalalignment','middle',...
    'horizontalalignment','right');
end
for expr=expr_min+0.5:1:expr_max-0.5
  y = ykn_start + (ykn_size / expr_range) * (expr - expr_min);
  text(label_x,y,sprintf('%d',expr),'fontsize',fontsize,'verticalalignment','middle',...
    'horizontalalignment','right');
end

% WRITE TITLES

if ~P.suppress_title
  title(upper(gene),'fontsize',fontsize+4);
end

text(0.5,label_y*2,'Copy number','fontsize',fontsize+2,'horizontalalignment','center');
text(label_x*3,0.5,P.expression_axis_label,'fontsize',fontsize+2,'rotation',90,...
  'horizontalalignment','center','verticalalignment','bottom');

% SET X and Y positions based on copy number and expression data

D.x = nan(slength(D),1);
D.y = D.x;
for i=1:slength(D)
  if D.have_cn(i)
    D.x(i) = xkn_start + xkn_size * (D.cn(i) - cn_min + (1-xscatter)/2 + xscatter*rand) / cn_ncats;
  else
    if P.show_unknowns
      D.x(i) = xunk_start + ((1-xscatter)/2 + xscatter*rand) * xunk_size;
    end
  end
  if D.have_expr(i)
    D.y(i) = ykn_start + ykn_size * (D.expr(i) - expr_min) / expr_range;
  else
    if P.show_unknowns
      D.y(i) = yunk_start + ((1-yscatter)/2 + yscatter*rand) * yunk_size;
    end
  end
end

% SUPERIMPOSE METHYLATION AND MUTATION INFO

% FACE COLOR = methylation
%    grey = no methylation data
%    blue->purple->red = low->high mutation (0->0.5->1)

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

% truncated_mut_color = yellow;
% other_mut_color = green;
% methylated_scheme = [1 0 0;0 0 0;0 0 1];

truncated_mut_color = red;
other_mut_color = pink;
methylated_scheme = [0.9 0.7 0;0 0.3 0;0 0 0.8];

D.facecolor = cell(slength(D),1);
D.edgecolor = cell(slength(D),1);
D.style = cell(slength(D),1);

for i = 1:slength(D)
  % face color = methylation
  if D.have_methyl(i)
    D.facecolor{i} = [D.methyl(i) 1 1-D.methyl(i)] * methylated_scheme;
  else  % no methylation data
    D.facecolor{i} = dkgrey;
  end
  % edge color = type of mutation
  if D.have_seq(i)
    if D.n_trunc(i), D.edgecolor{i} = truncated_mut_color;
    elseif D.n_other(i), D.edgecolor{i} = other_mut_color;
    else D.edgecolor{i} = ltgrey;
    end
    % size = number of mutations
    if D.n_tot(i)<2, D.style{i} = 'small';
    else D.style{i} = 'large';
    end
  else % not sequenced
    D.style{i} = 'bare';
    D.edgecolor{i} = D.facecolor{i};
  end
end

% SHAPE = cluster

D.shape = repmat({'cir'},slength(D),1);
D.shape(D.cluster==1) = repmat({'sqr'},sum(D.cluster==1),1);
D.shape(D.cluster==2) = repmat({'inv'},sum(D.cluster==2),1);
D.shape(D.cluster==3) = repmat({'tri'},sum(D.cluster==3),1);
D.shape(D.cluster==4) = repmat({'dmd'},sum(D.cluster==4),1);

% SORT POINTS so informative ones are displayed in front

D.methyl_score = nan(M.np,1);
idx = find(D.have_methyl);
if ~isempty(idx)
  D.methyl_score(idx) = (D.methyl(idx)-mean(D.methyl(idx))).^2;
end

D = sort_struct(D,{'n_tot','have_methyl','methyl_score'});

% DRAW DATA POINTS

for i=1:slength(D)
  if isnan(D.x(i))|isnan(D.y(i)), continue, end;
  draw_shape(D.shape{i},D.x(i),D.y(i),D.facecolor{i},D.style{i},D.edgecolor{i});
  if D.hypermutated(i)
    draw_shape('dot',D.x(i),D.y(i),red);
  elseif D.treated(i)
    draw_shape('dot',D.x(i),D.y(i),black);
  end
  if slength(D)-i < P.label_top_samples, text(D.x(i)-0.02,D.y(i)+0.02,D.name{i}(end-3:end),'fontsize',2); end
end

% CREATE IMAGEMAP

maptext=[];
for i=slength(D):-1:1
  if isnan(D.x(i))|isnan(D.y(i)), continue, end;
  if strcmp(D.style{i},'small'), d = 0.02;
  elseif strcmp(D.style{i},'large'), d = 0.03;
  else d = 0.016; end
  x = (axispos(1) + (D.x(i)*axispos(3))) * (paperpos(3)*paperres);
  y = (1-(axispos(2)+(D.y(i)*axispos(4)))) * (paperpos(4)*paperres);
  r = d/2 * (paperpos(3)*paperres);
  maptext = [maptext '<area title="' D.name{i} '" shape="circ" coords="' ...
    num2str(round(x)) ',' num2str(round(y)) ',' num2str(round(r)) '">' char(10)];
end

hold off
success = true;

