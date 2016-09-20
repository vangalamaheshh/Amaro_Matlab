function display_mutation_bargraph2(D,P)
% Mike Lawrence 2010

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'maxcount',[]);
P=impose_default_value(P,'maxcov',[]);
P=impose_default_value(P,'maxrate',[]);
P=impose_default_value(P,'barwidth',0.8);
if isfield(P,'off_scale_xadj'), P.text_xadj = P.off_scale_xadj; end
P=impose_default_value(P,'text_xadj',[]);
P=impose_default_value(P,'plot_width',0.95);
P=impose_default_value(P,'plot_left',0.03);
P=impose_default_value(P,'plot_height',0.28);
P=impose_default_value(P,'plot_ypos',[0.70 0.39 0.07]);
P=impose_default_value(P,'y_axis_fontsize',10);
P=impose_default_value(P,'ylabel_fontsize',13);
P=impose_default_value(P,'off_scale_fontsize',10);
P=impose_default_value(P,'sample_name_fontsize',10);

if isfield(D,'display_name')
  names = D.display_name;
elseif isfield(D,'names')
  names = D.names;
elseif isfield(D,'name')
  names = D.name;
else
  error('D contains neither name, nor names, nor display_name');
end

% try to shorten names
tmp = regexprep(names,'.*-(\d+)(-.*)?','$1');
if ~any(strcmp(tmp,names)) && length(unique(names))==length(unique(tmp)), names=tmp; end


% extract data
nnon_tot = D.nnon_tot;
nsil_tot = D.nsil_tot;
rate_sil = D.rate_sil;
rate_non = D.rate_non;
N_tot = D.N_tot;

if isfield(D,'ndbsnp_tot')
  ndbsnp_tot = D.ndbsnp_tot;
end
if isfield(D,'rate_dbsnp')
  rate_dbsnp = D.rate_dbsnp;
end

%%%%%%%  SET UP PLOT

set(gcf,'Color',[0.55 0.55 0.55],'InvertHardcopy','off')  % InvertHardcopy=off: keep grey bkgd in printed version


if ~exist('ndbsnp_tot','var')
  cmap = [0 0 0.5625;0.5 0 0];   % dk.blue = nonsilent; dk.red = silent;
else
  cmap = [0 0 0.5625;0.5 0 0;0.3 1 0.3];   % dk.blue = nonsilent; dk.red = silent; green = dbSNP
end
colormap(cmap);

left=P.plot_left; width=P.plot_width;
height=P.plot_height;
ypos=P.plot_ypos;
yadj = 0.9;
xmarg = 0.8;
if isempty(P.text_xadj)
  if length(names)<30, xadj = 0.1; elseif length(names)<100, xadj = 0.2; else xadj = 0.33; end
else
  xadj = P.text_xadj;
end
bar_params = {P.barwidth,'edgecolor','none'};
bar_stacked_params = {P.barwidth,'stacked','edgecolor','none'};
axis_params = {'fontsize',P.y_axis_fontsize};

%%%%   background colors to show sample type

if isfield(D,'type') & ~isfield(D,'label'), D = rename_field(D,'type','label'); end
if isfield(D,'label')
  D = make_numeric(D,'label');
  axes('position',[left 0 width 1],'visible','off');
  xlim([1-xmarg length(names)+xmarg]);
  labelcolors = distinct_colors(nan);
  for i=1:slength(D)
    if D.label(i)>=1 && D.label(i)<=size(labelcolors,1)
      col = labelcolors(D.label(i),:);
      pos = [i-0.4 0 1 1];    % maybe change i-0.4 to i-xadj
      rectangle('position',pos,'edgecolor',col,'facecolor',col);
end,end,end

%%%%%%%%%%%%  COUNTS

axes('position',[left ypos(1) width height]);

if ~exist('ndbsnp_tot','var')
  bar([nnon_tot nsil_tot],bar_stacked_params{:})
  leg = {'nonsilent','silent'};
  y = nnon_tot+nsil_tot;
else
  bar([nnon_tot nsil_tot ndbsnp_tot],bar_stacked_params{:})
  leg = {'nonsilent','silent','dbSNP'};
  y = nnon_tot+nsil_tot+ndbsnp_tot;
end

if isempty(P.maxcount)
  P.maxcount = min(1.1*nanmax(y),20*nanmedian(y));
end

xlim([1-xmarg length(names)+xmarg]);
set(gca,'xtick',[],axis_params{:});
if ~isempty(P.maxcount) && P.maxcount>0
  ylim([0 P.maxcount]);
  for i=1:length(names)
    if y(i)>P.maxcount
      text(i+xadj,P.maxcount*yadj,num2str(round(y(i))),'horizontalalignment','center','rotation',90,...
            'verticalalignment','middle','color',[1 1 1],'fontsize',P.off_scale_fontsize);
    end
  end
end
ylabel('# mutations','fontsize',P.ylabel_fontsize);

% MUTATION TYPE LEGEND
%legend('nonsilent','silent','dbSNP','location',P.legend);
xx = (length(names)+xmarg) * 0.86;
for i=1:length(leg)
  yy = P.maxcount * (1-0.1*(length(leg)-i+1));
  text(xx,yy,leg{i},'color',cmap(i,:),'fontweight','bold');
end  

%%%%%%%%%%%%  COVERAGE

axes('position',[left ypos(2) width height]);

y = N_tot/1e6;
bar(y,bar_params{:});
if isfield(D,'coverage_bar_colors')
  hold on
  q = nan(size(N_tot));
  for i=1:slength(D)
    qq = q; qq(i) = y(i);
    bar(qq,'edgecolor',D.coverage_bar_colors(i,:),'facecolor',D.coverage_bar_colors(i,:));
  end
  hold off
end

if isempty(P.maxcov)
  P.maxcov = min(1.1*nanmax(y),20*nanmedian(y));
end

xlim([1-xmarg length(names)+xmarg]);
set(gca,'xtick',[],axis_params{:});
if ~isempty(P.maxcov) && P.maxcov>0
  ylim([0 P.maxcov]);
  for i=1:length(names)
    if y(i)>P.maxcov
      text(i+xadj,P.maxcov*yadj,num2str(round(y(i))),'horizontalalignment','center','rotation',90,...
            'verticalalignment','middle','color',[1 1 1],'fontsize',P.off_scale_fontsize);
    end
  end
end
ylabel('Mb sequenced','fontsize',P.ylabel_fontsize);

% SAMPLE TYPE LEGEND

if isfield(D,'label') && isfield(P,'label_names') && exist('labelcolors','var')
  nl = length(P.label_names)
  tx = slength(D) * 0.03;
  ty = P.maxcov * 0.08;
  for i=nl:-1:1
    text(tx,ty,P.label_names{i},'BackgroundColor',labelcolors(i,:),'interpreter','none');
    ty=ty+(P.maxcov*0.18);
    if ty>P.maxcov*0.5
      tx = tx + (slength(D)*0.25);
      ty = P.maxcov * 0.08;
    end
  end
end

%%%%%%%%%%%%  RATES

axes('position',[left ypos(3) width height]);

if ~exist('rate_dbsnp','var')
  bar(1e6*[rate_non rate_sil],bar_stacked_params{:})
  y = 1e6*(rate_non+rate_sil);
else
  bar(1e6*[rate_non rate_sil rate_dbsnp],bar_stacked_params{:})
  y = 1e6*(rate_non+rate_sil+rate_dbsnp);
end

if isempty(P.maxrate)
  P.maxrate = min(1.1*nanmax(y),20*nanmedian(y));
end

xlim([1-xmarg length(names)+xmarg]);
if ~isempty(P.maxrate) && P.maxrate>0
  ylim([0 P.maxrate]);
  for i=1:length(names)
    if y(i)>P.maxrate
      text(i+xadj,P.maxrate*yadj,num2str(round(y(i))),'horizontalalignment','center','rotation',90,...
            'verticalalignment','middle','color',[1 1 1],'fontsize',P.off_scale_fontsize);
    end
  end
end
ylabel('mutations / Mb','fontsize',P.ylabel_fontsize);

xlim([0.2 length(names)+0.8]);
set(gca,'xtick',1:length(names),'xticklabel',names,axis_params{:});
set(gca,'TickLength',[0 0]);

xticklabel_rotate_simplified(xadj+(1:length(names)),90,names,'fontsize',P.sample_name_fontsize,'interpreter','none');
