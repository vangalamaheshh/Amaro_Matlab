function draw_mutation_rate_histogram(M,P)

% default parameters

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'exclude_these_genes',{});
P=impose_default_value(P,'mark_treated_and_secondary',false);
P=impose_default_value(P,'color_hypermutated_purple',false);
P=impose_default_value(P,'color_hypermutated_red',false);
P=impose_default_value(P,'label_hypermutated_samples',false);
P=impose_default_value(P,'manually_identify_hypermutated',{});
P=impose_default_value(P,'hist_with_error_bars',false);
P=impose_default_value(P,'hist_with_trend_lines',false);
P=impose_default_value(P,'filename','histogram.png');
P=impose_default_value(P,'axismax',200);
P=impose_default_value(P,'label_position_manual_adjustment',{});
P=impose_default_value(P,'output_table_of_rates',false);
P=impose_default_value(P,'table_of_rates_filename','sample_rates.txt');

P.color_hypermutated = P.color_hypermutated_purple | P.color_hypermutated_red;

fprintf('Drawing mutation-rate histogram\n');

% calculate total coverage for each sample
% and break down into silent and nonsilent "bases at risk"

hist_genes = 1:M.ng;
exclude_list = listmap(P.exclude_these_genes, M.gene.name);
hist_genes = setdiff(hist_genes, exclude_list);

N_breakdown = repmat(M.N_cov(:,M.TOT,:),1,24) .* repmat(M.breakdown.frac, [1 1 M.np]);
N_si = round(squeeze(sum(sum(N_breakdown(hist_genes,1:8,:),2))));
N_ns = round(squeeze(sum(sum(N_breakdown(hist_genes,9:24,:),2))));

n_ns = squeeze(sum(M.n_nonsilent(hist_genes,M.TOT,:),1));
n_si = squeeze(sum(M.n_silent(hist_genes,M.TOT,:),1));

% compute rates and confidence intervals

[phat_ns cdf_ns] = binofit(n_ns, N_ns, 0.05);
[phat_si cdf_si] = binofit(n_si, N_si, 0.05);

% take care of NaNs

idx = union(find(isnan(phat_ns)),find(isnan(phat_si)));
phat_ns(idx) = 0;
phat_si(idx) = 0;
cdf_ns(idx,:) = 0;
cdf_si(idx,:) = 0;

% output table of rates if desired

if P.output_table_of_rates
  out = fopen(P.table_of_rates_filename,'wt');
  fprintf(out,['sample\ttreated\thypermutated\tn_nonsilent\tN_nonsilent\trate_nonsilent' ...
               '\tn_silent\tN_silent\trate_silent\n']);
  for p=1:M.np
    fprintf(out, '%s\t%d\t%d\t%d\t%d\t%e\t%d\t%d\t%e\n', M.patient.name{p},...
       M.patient.treated(p), M.patient.hypermutated(p),...
       n_ns(p),N_ns(p), phat_ns(p), n_si(p), N_si(p), phat_si(p));
  end  
  fclose(out);
  fprintf('Type "return" to continue or "dbquit" to exit before outputting histogram\n');
  keyboard
end

if P.mark_treated_and_secondary
% OLD METHOD: NOW EVERYTHING IS LOADED DURING load_mutdata2()
%   SA = tab2struct('/xchip/tcga/gbm/results/sample_info/Sample_annotation_20080407.txt');

  treated = find(M.patient.treated);
  secondary = find(M.patient.secondary);
  sec_and_treated = intersect(secondary, treated);
  sec_or_treated = union(secondary, treated);
  sec_only = setdiff(secondary, sec_and_treated);
  treated_only = setdiff(treated, sec_and_treated);
  nonsec_nontreated = setdiff(1:M.np, sec_or_treated);
end

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

% put points in order so that specially colored points are out in front

if P.color_hypermutated || P.label_hypermutated_samples
  if ~isempty(P.manually_identify_hypermutated)
    hmut = P.manually_identify_hypermutated;
  else
    hmut = identify_hypermutated_samples(M,P);
  end
end

k=zeros(M.np,1);
for i=1:M.np
  if P.mark_treated_and_secondary
    if ismember(i,sec_and_treated), k(i)=2; end
    if ismember(i,treated_only), k(i)=1; end
  elseif P.color_hypermutated
    if ismember(M.patient.name{i},hmut), k(i)=1; end
  end
end

[tmp ord] = sort(k);

% plot data

for j=1:M.np
  i=ord(j);

  color = [0 0 0];

  if P.mark_treated_and_secondary
    if ismember(i,sec_and_treated), color = [1 0 0]; end
    if ismember(i,treated_only), color = [0 0 1]; end
    if ismember(i,sec_only), color = [0 1 0]; end
  end

  if P.color_hypermutated
    if ismember(M.patient.name{i},hmut)
      if P.color_hypermutated_red, color = [1 0 0];
      elseif P.color_hypermutated_purple, color = [1 0 1];
      end
    end
  end

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

  if P.label_hypermutated_samples
    pname = M.patient.name{i};
    if ismember(pname, hmut)
      xspc = 7*axismax/120;
      yspc = 3*axismax/120;
      x = xcen*sc-xspc; y=ycen*sc+yspc;
      for adj=1:size(P.label_position_manual_adjustment,1)
        if strcmp(pname,P.label_position_manual_adjustment{adj,1})
          x=x+P.label_position_manual_adjustment{adj,2}(1);
          y=y+P.label_position_manual_adjustment{adj,2}(2);
        end
      end
      text(x, y, pname);
    end
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

[dev res] = interpret_print_filename(P.filename);
fprintf('Outputting histogram to %s\n', P.filename);
print(['-d' dev], ['-r' num2str(res)], P.filename);

end





function old_style_histogram()

error('OLD-STYLE HISTOGRAM IS NO LONGER FUNCTIONAL');

%%%%%%%%%% CODE NEEDS TO BE COMPLETELY UPDATED

close all;
hold on
set(gca,'YDir','reverse');
ylabel('nonsilent mutations / MB-at-risk');
xlabel('silent mutations / MB-at-risk');
scatter(si(nonsec_nontreated),ns(nonsec_nontreated),'k.');
scatter(si(treated_only),ns(treated_only),'bx');
scatter(si(sec_only),ns(sec_only),'gx');
scatter(si(sec_and_treated),ns(sec_and_treated),'rx');

% draw 1:1 line

line([0 100], [0 100],'Color',[0 0 0]);

% fit to y = 1x + b

x=si;
y=ns;
a=1;
b = [-10:0.01:10];
ssd = zeros(length(b),1);
for i=1:length(b);
  y_model = (a*x)+b(i);
  sd = (y - y_model) .^ 2;
  ssd(i) = sum(sd);
end

b_best = b(find(ssd==min(ssd)));

% draw best-fit line
%   b=0.4 with all genes included
%   b=-0.9    with known cancer genes removed

line([0 100], [0 100]+b_best,'Color',[0 1 0]);

set(gca,'XLim',[0 100]);
set(gca,'YLim',[0 100]);

hold off
print_D(hist_png,{{'png','-r180'}});

return




function even_older_style_histogram()

error('EVEN-OLDER-STYLE HISTOGRAM IS NO LONGER FUNCTIONAL');

%%%%%%%%%% CODE NEEDS TO BE COMPLETELY UPDATED

hist_genes = 1:ng;
exclude_list = listmap(hist_excludes_these_genes, gene_name);
hist_genes = hist_genes(~ismember(genes_to_use,exclude_list));

pf =  squeeze(sum(n_potfun(hist_genes,TOT,:),1));
s = squeeze(sum(n_silent(hist_genes,TOT,:),1));
c = squeeze(sum(N_cov(hist_genes,TOT,:),1));

c = c./1000000;
pf = pf./c;
s = s./c;   % muts per MB

tot = pf;
[idx res] = outliers(tot,'iqr');

% find where MSH6 muts are

msh6 =  squeeze(n_potfun(find(strcmp(gene_name(genes_to_use),'MSH6')),TOT,:));
msh2 =  squeeze(n_potfun(find(strcmp(gene_name(genes_to_use),'MSH2')),TOT,:));

% find which samples were "treated"

SA = tab2struct('/xchip/tcga/gbm/results/sample_info/Sample_annotation_20080331.txt');
treated_list = SA.collaboratorParticipantId(grep('x',SA.treated,1));

% lists of secondary-tumor samples

sec_list = {;...
'TCGA-02-0010';...
'TCGA-02-0028';...
'TCGA-02-0102';...
'TCGA-02-0114';...
'TCGA-08-0525';...
};

% generate figure

tmp=[ round(pf)'; round(s)'];
xx=zeros(max(tmp(1,:))+3,max(tmp(2,:))+3);
x2=xx;
x6=xx;
xhm=xx;
xis150=xx;
xtreated=xx;
xoutlier=xx;
xis22=xx;
xissec=xx;
xistreated_and_sec = xx;
for i=1:size(tmp,2)
  a=tmp(1,i)+2;
  b=tmp(2,i)+2;
  xx(a,b)=xx(a,b)+1;
  x2(a,b)=x2(a,b)+msh2(i);
  x6(a,b)=x6(a,b)+msh6(i);
  xhm(a,b)=xhm(a,b)+ismember(patient_name{patients_to_use(i)},old_hypermutated_list);

  xis150(a,b)=xis150(a,b)+strcmp(patient_name{patients_to_use(i)},'TCGA-06-0150');

  treated=ismember(patient_name{patients_to_use(i)},treated_list);
  xtreated(a,b)=xtreated(a,b)+treated;

  xoutlier(a,b)=xoutlier(a,b)+ismember(i,idx);
  xis22(a,b)=xis22(a,b)+(i==22);

  issec=ismember(patient_name{patients_to_use(i)},sec_list);
  xissec(a,b)=xissec(a,b)+issec;

  if issec & treated
      xissec(a,b)=xissec(a,b)-1;
      xtreated(a,b)=xtreated(a,b)-1;
      xistreated_and_sec(a,b)=xistreated_and_sec(a,b)+1;
  end

end

close all
imagesc(xx)
colorbar

if 0
for y=1:size(xx,1)
  for x=1:size(xx,2)
    if xoutlier(y,x)
      if ~xis22(y,x) color=[0.4 0.4 0.4];
      else color=[0.2 0.2 0.2]; end
      rectangle('position',[x-.5 y-.5 1 1],'EdgeColor',color);
    end
    if xhm(y,x)
      text(x-0.5,y,'x  ','color',[1 0 0]);
    end
    if x2(y,x)
      text(x-0.3,y,repmat('x',1,x2(y,x)),'color',[0 0 0]);
    end
    if x6(y,x)
      text(x+0.16,y,repmat('x',1,x6(y,x)),'color',[1 1 1]);
    end
    if xis150(y,x)
      text(x-0.3,y,'150','color',[0 0 0]);
    end
    if xtreated(y,x)
      text(x-0.5,y,repmat('- ',1,xtreated(y,x)),'color',[1 .5 .5]);
    end
  end
end
end

for y=1:size(xx,1)
  for x=1:size(xx,2)
%    if xoutlier(y,x)
%      color=[0.4 0.4 0.4];
%      rectangle('position',[x-.5 y-.5 1 1],'EdgeColor',color);
%    end
    if xistreated_and_sec(y,x)
      text(x+0.3-0.1*(xistreated_and_sec(y,x)),y,repmat('x',1,xistreated_and_sec(y,x)),'color',[1 0 0]);
    elseif xtreated(y,x)
      text(x+0.15-0.1*(xtreated(y,x)),y,repmat('x',1,xtreated(y,x)),'color',[1 1 1]);
    elseif xissec(y,x)
      text(x-0.1*(xissec(y,x)),y,repmat('x',1,xissec(y,x)),'color',[0 0  1]);
    end
 end
end


ylabel('nonsilent mutations/MB');
xlabel('silent mutations/MB');
set(gca,'Xtick',1:size(xx,2),'XtickLabel',['  ';num2str((0:(size(xx,2)-3))');'  ']);
set(gca,'Ytick',1:size(xx,1),'YtickLabel',['  ';num2str((0:(size(xx,1)-3))');'  ']);
print_D(hist_png,{{'png','-r180'}});

end


end
