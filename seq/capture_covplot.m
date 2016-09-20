function [sidx tidx] = capture_covplot(C,P)
% Mike Lawrence 2009-2010

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'omit_these_samples',[]);
P=impose_default_value(P,'omit_these_targets',[]);
P=impose_default_value(P,'omit_these_genes',[]);
P=impose_default_value(P,'unit','target');
P=impose_default_value(P,'x_lastsort',[]);
P=impose_default_value(P,'y_lastsort',[]);
P=impose_default_value(P,'impute_full_coverage',false);

% SELECT SAMPLES
if ~isempty(P.omit_these_samples)
  keep = setdiff(1:C.ns,P.omit_these_samples);
  C.sample = reorder_struct(C.sample,keep);
  C.ns = length(keep);
%  C.cov = C.cov(:,keep);
  if isfield(C,'fcov'), C.fcov = C.fcov(:,keep); end
  if isfield(C,'gcov'), C.gcov = C.gcov(:,keep); end
  if isfield(C,'fgcov'), C.fgcov = C.fgcov(:,keep); end
  if ~isempty(P.x_lastsort), P.x_lastsort = P.x_lastsort(keep); end
  if isfield(P,'mutation_sample_names_to_match'),
    P.mutation_sample_names_to_match = P.mutation_sample_names_to_match(keep);
  end
end

% SELECT TARGETS
if ~isempty(P.omit_these_targets)
  keep = setdiff(1:C.nt,P.omit_these_targets);
  C.targ = reorder_struct(C.targ,keep);
  C.nt = length(keep);
%  C.cov = C.cov(keep,:);
  if isfield(C,'fcov'), C.fcov = C.fcov(keep,:); end
  if strcmp(P.unit,'target') && ~isempty(P.y_lastsort) , P.y_lastsort = P.y_lastsort(keep); end
end

% SELECT GENES
if ~isempty(P.omit_these_genes)
  if ~isempty(P.genelist), error('P.omit_these_genes and P.genelist cannot both be specified'); end
  keep = setdiff(1:C.ng,P.omit_these_genes);
  C.gene = reorder_struct(C.gene,keep);
  C.ng = length(keep);
  if isfield(C,'gcov'), C.gcov = C.gcov(keep,:); end
  if isfield(C,'fgcov'), C.fgcov = C.fgcov(keep,:); end
  if strcmp(P.unit,'gene')  && ~isempty(P.y_lastsort) , P.y_lastsort = P.y_lastsort(keep); end
end

P=impose_default_value(P,'genelist',[]);
P=impose_default_value(P,'genelist_fontsize',[]);
P=impose_default_value(P,'y_secondary_sort_by_gc',true);
P=impose_default_value(P,'dont_sort_these_samples',[]);
P=impose_default_value(P,'vertical_xlabels',true);
P=impose_default_value(P,'title','');
P=impose_default_value(P,'show_colorbar',true);
P=impose_default_value(P,'colorbar_fontsize',10);
P=impose_default_value(P,'omit_sample_names',slength(C.sample)>150);
if isfield(C.sample,'short'), nn= C.sample.short;
elseif isfield(C.sample,'name'), nn = C.sample.name;
else error('Can''t find sample names');
end
% try to shorten names
if iscellstr(nn)
  tmp = regexprep(nn,'.*-(\d+)(-.*)?','$1');
  if ~any(strcmp(tmp,nn)) && length(unique(nn))==length(unique(tmp)), nn=tmp; end
end

P=impose_default_value(P,'sample_names',nn);
P=impose_default_value(P,'sample_name_fontsize',10);
P=impose_default_value(P,'sample_name_xadj',0.33);
P=impose_default_value(P,'show_GC_plot',true);
P=impose_default_value(P,'GC_plot_fontsize',10);
P=impose_default_value(P,'mainheight',0.71);
P=impose_default_value(P,'mainypos',0.25);
P=impose_default_value(P,'mainwidth',0.94);
P=impose_default_value(P,'GC_subplot_width',0.05 + 0.05*(slength(C.sample)<30));
P=impose_default_value(P,'hide_axis_marks',true);
P=impose_default_value(P,'show_read_counts',...
     isfield(C.sample,'treads_ontarget') && ~all(isnan(C.sample.treads_ontarget)));
P=impose_default_value(P,'show_lane_counts',...
     isfield(C.sample,'tlanes') && ~all(isnan(C.sample.tlanes)));
P=impose_default_value(P,'lanes_axis_max',[]);
P=impose_default_value(P,'reads_axis_max',[]);
P=impose_default_value(P,'y_axis_lastsort_labels',[]);
P=impose_default_value(P,'y_axis_lastsort_labels_rotation',[]);
P=impose_default_value(P,'y_axis_lastsort_labels_fontsize',10);
P=impose_default_value(P,'x_axis_lastsort_labels',[]);
P=impose_default_value(P,'y_axis_lastsort_colors',distinct_colors(nan));
P=impose_default_value(P,'x_axis_lastsort_labels_fontsize',10);
P=impose_default_value(P,'x_axis_lastsort_colors',distinct_colors(nan));
if (~isempty(P.x_lastsort) && (length(P.x_axis_lastsort_colors) < max(P.x_lastsort) ||...
                               length(P.x_axis_lastsort_labels) < max(P.x_lastsort))) ||...
   (~isempty(P.y_lastsort) && (length(P.y_axis_lastsort_colors) < max(P.y_lastsort) ||...
                                length(P.y_axis_lastsort_labels) < max(P.y_lastsort))) ...
      error('Not enough x- or y-axis lastsort colors/labels');
end
P=impose_default_value(P,'bargraph_style','bidirectional');  % 'stacked' or 'bidirectional'
P=impose_default_value(P,'mutations',[]);
P=impose_default_value(P,'superimpose_mutations',true);
P=impose_default_value(P,'mutation_sample_names_to_match',P.sample_names);
P=impose_default_value(P,'mutations_match_margin',10);
if isfield(P,'dump_filename')
  fprintf('Please use P.dump_cov_filename instead of P.dump_filename\n');
  P = rename_field(P,'dump_filename','dump_cov_filename');
end
P=impose_default_value(P,'dump_cov_filename',[]);
P=impose_default_value(P,'dump_mut_filename',[]);

% SELECT TARGET or GENE data, copy to D

if strcmp(P.unit,'target')
  if isempty(P.genelist), idx = 1:C.nt;
  else
    idx=[];
    zz = listmap(C.targ.gene,P.genelist);
    for i=1:length(P.genelist)
      idx = [idx; find(zz==i)];
    end
    %%%% doesn't matter... alphabetical sort is still being imposed below (tidx3)
    %idx = find(ismember(C.targ.gene,P.genelist));
  end
  if isfield(C,'fcov_coding')
    D = C.fcov_coding(idx,:);
  elseif isfield(C,'fcov')
    D = C.fcov(idx,:);
  else
    D = [];
  end
  Y = reorder_struct(C.targ,idx);
elseif strcmp(P.unit,'gene')
  if isempty(P.genelist), idx = 1:C.ng;
  else
    idx=[];
    zz = listmap(C.gene.name,P.genelist);
    for i=1:length(P.genelist)
      idx = [idx; find(zz==i)];
    end
    %idx = find(ismember(C.gene.name,P.genelist));
  end
  if isfield(C,'fgcov_coding')
    D = C.fgcov_coding(idx,:);
  elseif isfield(C,'fgcov')
    D = C.fgcov(idx,:);
  else
    D = [];
  end
  Y = reorder_struct(C.gene,idx);
else
  error('P.unit must be "target" or "gene"');
end
ny = slength(Y);

clf,hold on
set(gcf,'Color',[0.55 0.55 0.55],'InvertHardcopy','off')  % InvertHardcopy=off: keep grey bkgd in printed version
%set(gcf,'position',[638    47   708   774]);

if ~P.impute_full_coverage
  % SORT SAMPLES
  if ~isempty(D)
    [tmp sidx1] = sort(mean(D,1),'descend');
  else
    sidx1 = 1:C.ns;
  end
  if ~isempty(P.x_lastsort)
    [tmp sidx2] = sort(P.x_lastsort(sidx1));
  else
    sidx2 = 1:C.ns;
  end
  sidx = sidx1(sidx2);
else
  sidx = 1:C.ns;
end

% SORT GENES/TARGETS

ss = setdiff(1:C.ns,P.dont_sort_these_samples);

if P.y_secondary_sort_by_gc && ~isempty(D)
  q = Y.gc;
  z = mean(D(:,ss),2);
  q(z>0)=nan;
  [tmp tidx1] = sort(q,'ascend');
else
  tidx1 = 1:ny;
end
if ~isempty(D)
  [tmp tidx2] = sort(mean(D(tidx1,ss),2),'descend');
else
  tidx2 = 1:ny;
end
if strcmp(P.unit,'target') && ~isempty(P.genelist)
  [tmp tidx3] = sort(Y.gene(tidx1(tidx2)));
else
  tidx3 = 1:ny;
end
if ~isempty(P.y_lastsort)
  [tmp tidx4] = sort(P.y_lastsort(tidx1(tidx2(tidx3))));
else
  tidx4 = 1:ny;
end
tidx = tidx1(tidx2(tidx3(tidx4)));

% MAIN COVERAGE PLOT
q = ~P.show_lane_counts + ~P.show_read_counts;
mainheight = P.mainheight+0.09*q; mainypos = P.mainypos-0.09*q;
if ~P.show_GC_plot, P.GC_subplot_width=0; mainxpos = 0.04;
else mainxpos = P.GC_subplot_width + 0.05; end
mainwidth = P.mainwidth - P.GC_subplot_width;
subplot('position',[mainxpos mainypos mainwidth mainheight]);
if ~isempty(D)
  imagesc(D(tidx,sidx),[0 1])
else
  fprintf('capture_covplot: imputing full coverage.\n');
  imagesc(1,[0 1]); 
  rectangle('position',[0.5 0.5 length(sidx) length(tidx)],'edgecolor',[0 0 0],'facecolor',[1 1 1]);
  xlim([0.5 length(sidx)+0.5]); ylim([0.5 length(tidx)+0.5]);
  text(0.5+length(sidx)/2,0.5+length(tidx)/2,'(imputed full coverage)','horizontalalign','center');
end
colormap(hot)
if P.show_colorbar
  params = {};
  if ~isempty(P.colorbar_fontsize), params = [params 'fontsize' P.colorbar_fontsize]; end
  colorbar(params{:})
end
title(P.title,'fontsize',20);

% X-AXIS LABELS
if ~isempty(P.x_axis_lastsort_labels)
  if isempty(P.x_lastsort)
    error('Must specify x_lastsort if specifying x_axis_lastsort_labels');
  else
    for i=1:length(P.x_axis_lastsort_labels)
      idx = find(P.x_lastsort(sidx)==i);
      if isempty(idx), continue; end
      xtick = median(idx);
      xmin = min(idx);
      xmax = max(idx);
      if ny<1000
        yyy = [1.01 -0.005];
        yyy2 = 1.022;
      else
        yyy = [1.002 -0.0115];
        yyy2 = 1.0135;
      end
      for yy=yyy
        rectangle('position',[xmin-0.5 yy*ny xmax-xmin+1 0.012*ny],...
          'facecolor',P.x_axis_lastsort_colors(i,:),...
          'edgecolor',P.x_axis_lastsort_colors(i,:),...
          'clipping','off');
      end
      text(xtick,yyy2*ny,P.x_axis_lastsort_labels{i},'fontsize',P.x_axis_lastsort_labels_fontsize,...
        'horizontalalign','center','interpreter','none');
      if i>1, line((xmin-0.5)*[1 1],ny*[yyy(2) yyy(1)+0.012],'color',[0 0 0],'clipping','off'); end
    end
  end
end

if P.omit_sample_names || P.impute_full_coverage
  set (gca,'xtick',[]);
  xlabel('sample');
%  text(C.ns/2,ny*1.03,'sample','fontsize',20);
else
  fsz = P.sample_name_fontsize;
  xadj = P.sample_name_xadj;  % 0.33
  if P.vertical_xlabels && C.ns>1
    xticklabel_rotate(xadj+[1:C.ns],90,P.sample_names(sidx),'fontsize',fsz,'interpreter','none');
  else
    set(gca,'xtick',xadj+[1:C.ns],'xticklabel',P.sample_names(sidx),'fontsize',fsz);
  end
end

% Y-AXIS LABELS
if ~isempty(P.genelist)
  ytick = [];
  if strcmp(P.unit,'gene'), g = Y.name(tidx);
  else g = Y.gene(tidx); end
  if length(P.genelist)<=6, rot=90;txsz=12; else rot=0;txsz=8; end
  if ~isempty(P.genelist_fontsize), txsz = P.genelist_fontsize; end
  for i=1:length(P.genelist)
    idx = find(strcmp(g,P.genelist{i}));
    thisy = (min(idx)+max(idx))/2;
    thisy = thisy + 0.025*ny/(length(idx));
    ytick = [ytick;thisy];
    if C.ns<100, ccx = 0.015; else ccx = 0.025; end 
    text(-ccx*C.ns,thisy,P.genelist{i},'fontsize',txsz,'verticalalign','middle',...
      'horizontalalign','center','rotation',rot,'interpreter','none');
  end
%  [tmp ord] = sort(ytick);
%  set(gca,'ytick',ytick(ord),'yticklabel',P.genelist(ord),'fontsize',10);
elseif ~isempty(P.y_axis_lastsort_labels)
  if isempty(P.y_lastsort)
    error('Must specify y_lastsort if specifying y_axis_lastsort_labels');
  else
    ytick = [];
    for i=1:length(P.y_axis_lastsort_labels)
      idx = find(P.y_lastsort(tidx)==i);
      if isempty(idx), continue; end
      ytick = [ytick;median(idx)];
      ymin = min(idx);
      ymax = max(idx);
      if C.ns<30
        xxx = [0.012 1.026];
        xxx2 = -0.01;
      elseif C.ns<80
        xxx = [0 1.012];
        xxx2 = -0.015;
      else
        xxx = [-0.008 1.008];
        xxx2 = -0.025;
      end
      for xx=xxx
        rectangle('position',[xx*C.ns ymin-0.5 0.01*C.ns ymax-ymin+1],...
          'facecolor',P.y_axis_lastsort_colors(i,:),...
          'edgecolor',P.y_axis_lastsort_colors(i,:),...
          'clipping','off');
      end
      if ~isempty(P.y_axis_lastsort_labels_rotation)
        if length(P.y_axis_lastsort_labels_rotation)>1, iii=i; else iii=1; end
        rot = P.y_axis_lastsort_labels_rotation(iii);
      else
        rot = 0;
      end
      text(xxx2*C.ns,median(idx),P.y_axis_lastsort_labels{i},'fontsize',...
           P.y_axis_lastsort_labels_fontsize,'rotation',rot,...
           'verticalalign','middle','horizontalalign','center','interpreter','none');
      if i>1, line(C.ns*[xxx(1) xxx(2)+0.01],(ymin-0.5)*[1 1],'color',[0 0 0],'clipping','off'); end
    end
%    [tmp ord] = sort(ytick);
%    set(gca,'ytick',ytick(ord),'yticklabel',P.y_axis_lastsort_labels(ord),'fontsize',14);
     set(gca,'ytick',[]);
  end 
else
%  text(-0.025*C.ns,ny/2,P.unit,'fontsize',20','verticalalign','middle','horizontalalign','center',...
%    'rotation',90,'interpreter','none');
  %ylabel(P.unit,'fontsize',20);

end

if P.hide_axis_marks
  set(gca,'visible','off');
end

mainplot_pos = get(gca,'position');

% process MUTATIONS if available, and (if specified), SUPERIMPOSE them on the graph

if ~isempty(P.mutations) && ~isempty(D)
  fprintf('Superimposing mutations\n');
  Dmut = zeros(size(D));
  M = P.mutations;
  require_fields(M,{'patient','chr','start','end'});
  M = make_numeric(M,{'start','end'});
  if ~isnumeric(M.chr), M.chr = convert_chr(M.chr); end
  for i=1:slength(M)
    idxa = find(strcmp(M.patient{i},P.mutation_sample_names_to_match),1);
    if isempty(idxa), continue; end
    idxb = find(Y.chr==M.chr(i) & Y.start<=M.end(i)+P.mutations_match_margin &...
                Y.end>=M.start(i)-P.mutations_match_margin,1);
    if isempty(idxb), continue; end
    Dmut(idxb,idxa) = Dmut(idxb,idxa) + 1;
    if P.superimpose_mutations
      xxx = find(sidx==idxa);
      yyy = find(tidx==idxb);
      pp = {'horizontalalignment','center','verticalalignment','middle'};
      if isfield(M,'color')  
        dkcolor = M.color(i,:); ltcolor = min(dkcolor+0.5,1);
      else
        ltcolor = [0.5 0.5 1];  dkcolor = [0 0 1];
      end
      text(xxx+0.12,yyy+0.12,'o','color',ltcolor,'fontsize',8,pp{:});
      text(xxx+0.12,yyy+0.12,'o','color',dkcolor,'fontsize',6,pp{:});
    end
  end    
end

% draw lines to divide genes (if genelist specified)
if ~isempty(P.genelist)
  tmp = find(~strcmp(Y.gene(tidx(2:end)),Y.gene(tidx(1:end-1))));
  for i=1:length(tmp)
    if C.ns<100, ccx=2; else ccx=5; end
    line([-ccx C.ns+ccx],[tmp(i) tmp(i)]+0.5,'color',[0 0 0],'clipping','off');
  end
end

% PLOT GC CONTENT
if P.show_GC_plot
  params = {};
  if ~isempty(P.GC_plot_fontsize), params = [params 'fontsize' P.GC_plot_fontsize]; end
  subplot('position',[0.014 mainplot_pos(2) P.GC_subplot_width mainplot_pos(4)]);
  if ny>1000, smooth_window = round(ny/1000); else smooth_window = 1; end
  if ny>1, yy1=1; yy2=ny; else yy1=0; yy2=1; end
  plot(smooth(Y.gc(flipud(tidx)),smooth_window),yy1:yy2); ylim([yy1 yy2]);
  xlabel('GC content',params{:});xlim([0.2 0.8]);
  set(gca,'ytick',[],'xtick',[0.2 0.5 0.8],params{:});
  line([0.5 0.5],[1 ny],'color',[0 0 0]);
end

if P.omit_sample_names
  subploty = 0.1;
else
  subploty = 0.03;
end
legendy = subploty-0.01;
grey = 0.6*ones(1,3);

if P.show_read_counts

% BARGRAPH on-target reads

subplot('position',[mainplot_pos(1) subploty mainplot_pos(3) 0.07]);
subploty = subploty + 0.03;
legendy = legendy - 0.05;
t = C.sample.treads_ontarget / 1e6;
n = C.sample.nreads_ontarget / 1e6;
if strcmp(P.bargraph_style,'stacked')
  bar([t(sidx),n(sidx)],0.6,'stacked');
  maxy = max(t+n)*1.1; if isnan(maxy) maxy = 1; end
  if ~isempty(P.reads_axis_max), maxy = P.reads_axis_max; end
  ylim([0 maxy]);
elseif strcmp(P.bargraph_style,'bidirectional')
  hold on
  bar(t(sidx),0.6,'facecolor',[0 0 0]);
  bar(-n(sidx),0.6,'facecolor',grey,'edgecolor',grey);
  maxy = max([t;n])*1.1; if isnan(maxy) maxy = 1; end
  if ~isempty(P.reads_axis_max), maxy = P.reads_axis_max; end
  ylim([-maxy maxy]);
  set(gca,'yticklabel',regexprep(char2cell(get(gca,'yticklabel')),'-|\s',''));
  for i=1:C.ns
    if t(sidx(i))>maxy, text(i+0.25,maxy/2,num2str(t(sidx(i))),'color',[1 1 1],'fontsize',9,...
      'rotation',90,'verticalalign','middle','horizontalalign','center'); end
    if n(sidx(i))>maxy, text(i+0.25,-maxy/2,num2str(n(sidx(i))),'color',[0 0 0],'fontsize',9,...
      'rotation',90,'verticalalign','middle','horizontalalign','center'); end
  end
else error('P.bargraph_style should be either "stacked" or "bidirectional"');
end
xlim(0.5+[0 C.ns]);
set(gca,'xtick',[],'xcolor',[1 1 1]);
ylabel({'on-target','reads','(millions)'});

end

if P.show_lane_counts

% BARGRAPH lanes sequenced
subplot('position',[mainplot_pos(1) subploty mainplot_pos(3) 0.07]);
t = C.sample.tlanes;
n = C.sample.nlanes;

if strcmp(P.bargraph_style,'stacked')
  bar([t(sidx),n(sidx)],0.6,'stacked');
  maxy = max(t+n)*1.1; if isnan(maxy) maxy = 1; end
  if ~isempty(P.lanes_axis_max), maxy = P.lanes_axis_max; end
  ylim([0 maxy]);
elseif strcmp(P.bargraph_style,'bidirectional')
  hold on
  bar(t(sidx),0.6,'facecolor',[0 0 0]);
  bar(-n(sidx),0.6,'facecolor',grey,'edgecolor',grey);
  maxy = max([t;n])*1.1; if isnan(maxy) maxy = 1; end
  if ~isempty(P.lanes_axis_max), maxy = P.lanes_axis_max; end
  ylim([-maxy maxy]);
  set(gca,'yticklabel',regexprep(char2cell(get(gca,'yticklabel')),'-|\s',''));
  for i=1:C.ns
    if t(sidx(i))>maxy, text(i+0.25,maxy/2,num2str(t(sidx(i))),'color',[1 1 1],'fontsize',9,...
      'rotation',90,'verticalalign','middle','horizontalalign','center'); end
    if n(sidx(i))>maxy, text(i+0.25,-maxy/2,num2str(n(sidx(i))),'color',[0 0 0],'fontsize',9,...
      'rotation',90,'verticalalign','middle','horizontalalign','center'); end
  end
else error('P.bargraph_style should be either "stacked" or "bidirectional"');
end
xlim(0.5+[0 C.ns]);
set(gca,'xtick',[],'xcolor',[1 1 1]);
ylabel('lanes');

end% ly = 0.04 if both; 0.07 if not

if P.show_lane_counts || P.show_read_counts
  % T/N legend
  subplot('position',[0.02 legendy 0.03 0.1]);
  set(gca,'visible','off');
  rectangle('position',[0 0 1 1],'facecolor',[1 1 1],'edgecolor',[0 0 0],'clipping','off');
  rectangle('position',[0.3 0.45 0.4 0.45],'facecolor',[0 0 0]);
  text(0.5,0.7,'T','fontsize',10,'color',[1 1 1],'horizontalalignment','center','verticalalignment','middle');
  line([0 1],[0.45 0.45],'color',[0 0 0]);
  rectangle('position',[0.3 0.15 0.4 0.3],'facecolor',grey,'edgecolor',grey);
  text(0.5,0.2,'N','fontsize',10,'color',[0 0 0],'horizontalalignment','center','verticalalignment','middle');
end
  
% DONE
hold off

% DUMP COVERAGE and/or MUTCOUNTS (if requested)

if ~isempty(P.dump_mut_filename) && ~isempty(D)
  Tmut = reorder_struct(Y,tidx);
  if isfield(C.sample,'filename'), nn = C.sample.filename;
  elseif isfield(C.sample,'med'), nn = C.sample.med;
  elseif isfield(C.sample,'name'), nn = C.sample.name;
  elseif isfield(C.sample,'short'), nn = C.sample.short;
  else error('Can''t find sample names');
  end
  fieldname = genfieldname(regexprep(nn,'-','_'));
  for i=1:C.ns
    Tmut = setfield(Tmut,fieldname{sidx(i)},Dmut(tidx,sidx(i)));
  end
  fprintf('Writing mutation counts to %s\n',P.dump_mut_filename);
  save_struct(Tmut,P.dump_mut_filename);
end

if ~isempty(P.dump_cov_filename) && ~isempty(D)
  T = reorder_struct(Y,tidx);
  if isfield(C.sample,'filename'), nn = C.sample.filename;
  elseif isfield(C.sample,'med'), nn = C.sample.med;
  elseif isfield(C.sample,'name'), nn = C.sample.name;
  elseif isfield(C.sample,'short'), nn = C.sample.short;
  else error('Can''t find sample names');
  end
  fieldname = genfieldname(regexprep(nn,'-','_'));
  for i=1:C.ns
    T = setfield(T,fieldname{sidx(i)},D(tidx,sidx(i)));
  end
  fprintf('Writing coverage to %s\n',P.dump_cov_filename);
  save_struct(T,P.dump_cov_filename);
end
