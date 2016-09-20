function draw_multichr_plot(S,seg,D,chrs,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'draw_grey_rects',false);
P = impose_default_value(P,'segfile',[]);
P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'left_clip',[]);
P = impose_default_value(P,'right_clip',[]);
P = impose_default_value(P,'inset_horiz_margin',0.06);
P = impose_default_value(P,'inset_chr',[]);
P = impose_default_value(P,'inset_left',[]);
P = impose_default_value(P,'inset_right',[]);
P = impose_default_value(P,'show_all_rearrangements_despite_inset',false);
P = impose_default_value(P,'show_inset_copy_number',false);
P = impose_default_value(P,'refseq',[]);
P = impose_default_value(P,'name_only_these_genes',[]);
P = impose_default_value(P,'gene_height',0.05);
P = impose_default_value(P,'midstrip',0.1+0.3*~isempty(P.inset_chr));
P = impose_default_value(P,'inset_gene_whichmidstrip',1);
P = impose_default_value(P,'inset_maxcnratio_dashedlinewidth',0.5);
P = impose_default_value(P,'show_segseq_dots',true);
P = impose_default_value(P,'show_segseq_segs',false);
P = impose_default_value(P,'line_anchors','segfile');
P = impose_default_value(P,'line_curvedness',0.1);
P = impose_default_value(P,'scatter_point_size',20);
P = impose_default_value(P,'dRanger_linewidth',2);
P = impose_default_value(P,'expanded_color_scheme',false);
%P = impose_default_value(P,'centromere_x_adjust',[]);
%P = impose_default_value(P,'centromere_y_adjust',[]);
P = impose_default_value(P,'centromere_height',0.035);
P = impose_default_value(P,'xon',zeros(length(chrs),1));
P = impose_default_value(P,'compression_factor',[]);
P = impose_default_value(P,'inset_compression_factor',P.compression_factor);
P = impose_default_value(P,'fontsize',16);
P = impose_default_value(P,'gene_name_fontsize',P.fontsize);
P = impose_default_value(P,'no_labeling',false);
P = impose_default_value(P,'note_max_copy_ratio',true);

nchrs = length(chrs);
if nchrs<1, error('Need at least one chr'); end
if length(unique(chrs))~=nchrs, error('duplicate chrs not allowed'); end

if length(P.xon)~=nchrs, error('xon must have one entry per chromosome if specified'); end

if ~isempty(S)
  sn = cell(nchrs,1);
  for c=1:nchrs
    sn{c} = reorder_struct(S.RATIOS,S.RATIOS.chr==chrs(c));
end,end
if ~isempty(seg)
  segn = cell(nchrs,1);
  for c=1:nchrs
    segn{c} = reorder_struct(seg.SEG,seg.SEG.chr==chrs(c));
end,end
if ~isempty(P.segfile)
  S = load_cnseg_file(P.segfile);
  S.copyratio = 2.^S.segmean;
  min_nprobes = 8;
  S = reorder_struct(S,S.nprobes>=min_nprobes);
  Sn = cell(nchrs,1);
  for c=1:nchrs
    Sn{c} = reorder_struct(S,S.chr==chrs(c));
end,end

% load centromere locations
try
  cen = load_cen(P.build);
catch me
  cen = [];
  fprintf('Failed to load centromere locations\n');
end

% load Refseq entry for inset gene (if specified)
if ~isempty(P.inset_chr)
  if length(P.inset_chr)>1, error('multiple inset chromosomes not yet supported'); end
  igchr = P.inset_chr;
  igc = find(chrs==igchr);
  if isempty(igc), error('chr%d is not on among specified chrs',P.inset_chr); end
  if isempty(P.refseq)
    R = load_refseq(P.build);
  else
    R = P.refseq;
  end
  R.chr = convert_chr(R.chr);
  Rc = reorder_struct(R,R.chr==P.inset_chr & R.tx_start<=P.inset_right & R.tx_end >=P.inset_left);
end

% determine copy-number dynamic range
v = cell(nchrs,1);
for c=1:nchrs
  v{c} = [];
  if exist('sn','var'), v{c} = [v{c}; sn{c}.ratios]; end
  if exist('Sn','var'), v{c} = [v{c}; Sn{c}.ratio]; end
end
v = cat(1,v{:});
mn = min(v);
mx = max(v);
md = (mn+mx)/2;

% determine chromosomal coordinates range
leftn = nan(nchrs,1);
rightn = nan(nchrs,1);
for c=1:nchrs
  pn = [];
  if exist('sn','var'), pn = [pn;double(sn{c}.windows)]; end
  if exist('segn','var'), pn = [pn;segn{c}.left;segn{c}.right]; end
  if exist('Sn','var'), pn = [pn;Sn{c}.start;Sn{c}.end]; end
  leftn(c) = min(pn);
  rightn(c) = max(pn);
end
% override?
if ~isempty(P.left_clip), leftn(~isnan(P.left_clip)) = P.left_clip(~isnan(P.left_clip)); end
if ~isempty(P.right_clip), rightn(~isnan(P.right_clip)) = P.right_clip(~isnan(P.right_clip)); end

left = min(leftn);
right = max(rightn);
width = right-left+1;

% SET UP DRAWING AREAS

% y areas
ymin = 0; ymax = 1;
topstrip = 0.05;
btmstrip = 0.05;
midstrip = P.midstrip;
chrstrip = (1-(topstrip+btmstrip+((nchrs-1)*midstrip)))/nchrs;
if (chrstrip<=0) error('margin strips too big'); end
if (mx-mn)>=2, compression_factor=1.2; else compression_factor = 1.2/(mx-mn); end
if ~isempty(P.compression_factor), compression_factor = P.compression_factor; end
ys = (0.5*chrstrip) / (compression_factor*(md-mn));
yn = nan(nchrs,1);
ytop = nan(nchrs,1);
y = 1-topstrip;
for c=1:nchrs
  ytop(c) = y;
  yn(c) = y - 0.5*chrstrip - (md*ys);
  y = y - chrstrip - midstrip;
end

% x-offsets
xon = P.xon;
label_x = double(width*0.01*1e-6);

% PLOT

clf;hold on
l2dir = P.line_curvedness;   % spline-flip toggle for line2()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INSET

if ~isempty(P.inset_chr)
  % set up gene inset
  if nchrs>2, error('need to modify code to take into account P.inset_gene_whichmidstrip'); end

  inset_horiz_margin = P.inset_horiz_margin;
  inset_left = 1e-6*(left+width*inset_horiz_margin);
  inset_right = 1e-6*(right-width*inset_horiz_margin);
  inset_width = inset_right - inset_left;

  inset_start = P.inset_left;
  inset_end = P.inset_right;
  inset_len = inset_end - inset_start + 1;
  isc = inset_width / inset_len;

  gene_height = P.gene_height;
  gene_mid = btmstrip+chrstrip+(0.5*midstrip);
  gene_top = gene_mid+(0.5*gene_height);
  gene_btm = gene_mid-(0.5*gene_height);

  ysi = ys / (P.inset_compression_factor / P.compression_factor);

  % draw "zoom lines" relating gene inset to its chromosome
  gx1 = 1e-6*(xon(igc)+inset_start); gx2 = 1e-6*(xon(igc)+inset_end);
  gy1 = yn(igc)+ysi*get_ratio(igc,inset_start); gy2 = yn(igc)+ysi*get_ratio(igc,inset_end);
  line([inset_left gx1],[gene_mid gy1],'color',[0.9 0.9 0.9],'linewidth',10);
  line([inset_right gx2],[gene_mid gy2],'color',[0.9 0.9 0.9],'linewidth',10);

 % draw intron/exon diagram
  line([inset_left inset_right],gene_mid*[1 1],'color',[0.5 0.5 0.5],'linewidth',4);
  already_named = {};
  for r=1:slength(Rc)
    rectangle('position',[inset_left+isc*(Rc.tx_start(r)-inset_start) gene_btm (isc*(Rc.tx_end(r)-Rc.tx_start(r))) gene_height],...
              'facecolor',[0.8 0.8 0.8]);
    if ~ismember(Rc.gene{r},already_named) && (isempty(P.name_only_these_genes) || ismember(Rc.gene{r},P.name_only_these_genes))
      text(inset_left+isc*(((Rc.tx_start(r)+Rc.tx_end(r))/2)-inset_start),gene_top,Rc.gene{r},...
           'horizontalalignment','center','fontsize',P.gene_name_fontsize,'verticalalignment','bottom');  % show gene name
      already_named{end+1} = Rc.gene{r};
    end
    for i=1:Rc.n_exons(r)
      st = Rc.exon_starts{r}(i); en = Rc.exon_ends{r}(i);
      rectangle('position',[inset_left+isc*(st-inset_start) gene_btm isc*(en-st) gene_height],'facecolor',[0 0 0]);
    end
  end

  % draw centromere
  if ~isempty(cen)
    cencolor = [0.6 0.6 1];
    cenh = 0.035;
    cleft = cen(chrs(igc),1); cright = cen(chrs(igc),2);
    if (cleft<inset_end && cright>inset_start)
      cratio = get_ratio(igc,(cleft+cright)/2);
      rectangle('position',[inset_left+isc*(cleft-inset_start) gene_btm isc*(cright-cleft) gene_height],...
                'curvature',1,'facecolor',cencolor,'edgecolor',cencolor);
    end
  end

  if P.show_inset_copy_number
    idx = find(sn{igc}.windows>=inset_start & sn{igc}.windows<=inset_end);
    scatter(inset_left+isc*(sn{igc}.windows(idx)-inset_start),...
            gene_top+2*gene_height+ysi*sn{igc}.ratios(idx),P.scatter_point_size,[0 0 0],'filled')
    % note max copy ratio
    [mxcnratio mxidx] = max(sn{igc}.ratios(idx));
    mxpos = sn{igc}.windows(idx(mxidx));
    dashedlinewidth = P.inset_maxcnratio_dashedlinewidth*inset_len;
    mxleft = mxpos - dashedlinewidth/2;
    mxright = mxpos + dashedlinewidth/2;
    line(inset_left+isc*([mxleft mxright]-inset_start),...
         (gene_top+2*gene_height+ysi*mxcnratio)*[1 1],'color',[0.5 0.5 0.5],'linestyle',':');
    text(inset_left+isc*(mxright-inset_start),...
         (gene_top+2*gene_height+ysi*mxcnratio),sprintf('CN ratio %-4.1f',mxcnratio),'interpreter','none');
  end
end


% draw grey rectangles (no longer supported)
if P.draw_grey_rects
  error('Need to fix this section to work with arbitrary number of chromosomes!');
  grey = ones(1,3)*0.8;
  rectangle('position',[left/1e6 1-topstrip width/1e6 topstrip],'facecolor',grey,'edgecolor',grey);
  rectangle('position',[left/1e6 btmstrip+chrstrip width/1e6 midstrip],'facecolor',grey,'edgecolor',grey);
  rectangle('position',[left/1e6 0 width/1e6 btmstrip],'facecolor',grey,'edgecolor',grey);
  
  if righta<right
    rectangle('position',[righta/1e6 0.5 (right-righta)/1e6 0.5],'facecolor',grey,'edgecolor',grey);
  end
  if rightb<right
    rectangle('position',[rightb/1e6 0 (right-rightb)/1e6 0.5],'facecolor',grey,'edgecolor',grey);
  end
end

% draw centromeres
if ~isempty(cen)
  cencolor = [0.6 0.6 1];
  cenh = P.centromere_height; % 0.035;
  for c=1:nchrs
    cleft = cen(chrs(c),1); cright = cen(chrs(c),2);
    cratio = get_ratio(c,(cleft+cright)/2);
    rectangle('position',[1e-6*(xon(c)+cleft) yn(c)+(ys*cratio-cenh/2) 1e-6*(cright-cleft) cenh],...
              'curvature',1,'facecolor',cencolor,'edgecolor',cencolor);
end,end


% draw segfile segments (if available)
if ~isempty(P.segfile) && exist('Sn','var')
  for c=1:nchrs
    for i=1:slength(Sn{c})
      line(1e-6*(xon(c)+[Sn{c}.start(i) Sn{c}.end(i)]),yn(c)+ys*Sn{c}.copyratio(i)*[1 1],'color',[0.6 0.6 0.6],'linewidth',4);
end,end,end

% segseq
if P.show_segseq_dots && exist('sn','var')
  for c=1:nchrs
    scatter(1e-6*(xon(c)+sn{c}.windows),yn(c)+ys*sn{c}.ratios,P.scatter_point_size,[0 0 0],'filled')
    if P.note_max_copy_ratio
      % note max copy ratio
      [mxcnratio mxidx] = max(sn{c}.ratios);
      mxpos = sn{c}.windows(mxidx);
      dashedlinewidth = width * 0.1;  %% line to extend over 10% of display width
      mxleft = mxpos - dashedlinewidth/2;
      mxright = mxpos + dashedlinewidth/2;
      line(1e-6*(xon(c)+[mxleft mxright]),(yn(c)+ys*mxcnratio)*[1 1],'color',[0.5 0.5 0.5],'linestyle',':');
      text(1e-6*(xon(c)+mxright),yn(c)+ys*mxcnratio,sprintf('CN ratio %-4.1f',mxcnratio),'interpreter','none');
    end
  end
end

if P.show_segseq_segs && exist('segn','var')
  for c=1:nchrs
    for i=1:slength(segn{c})
      line(1e-6*(xon(c)+[segn{c}.left(i),segn{c}.right(i)]),[1 1]*(yn(c)+ys*segn{c}.ratios(i)),'color','green','linewidth',2.5);
end,end,end

% dRanger rearrangements

% dRanger class colors
intercolor = [0.7 0 0.8];
intracolor = [0.1 0.7 0.1];
classes = {'inter_chr','deletion','inversion','long_range','tandem_dup'};
if P.expanded_color_scheme
  class_colors = {intercolor,[0 0 1],[0 1 1],[0.2 0.9 0.2],[1 0 0]};
else
  class_colors = {intercolor,intracolor,intracolor,intracolor,intracolor};
end
D.color = map_across(D.class,classes,class_colors);
drlw = P.dRanger_linewidth;   % dRanger linewidth

[tmp ord] = sort(listmap(D.class,classes));
D = reorder_struct(D,ord);
D.chr1 = convert_chr(D.chr1); D.chr2 = convert_chr(D.chr2);
D = make_numeric(D,{'pos1','pos2'});

if isempty(P.inset_chr) || P.show_all_rearrangements_despite_inset
  % if no inset gene, then show all rearrangements involving the two chromosomes
  for i=1:slength(D)
    c1 = find(D.chr1(i)==chrs);
    c2 = find(D.chr2(i)==chrs);
    if isempty(c1) || isempty(c2), continue; end

    p1 = D.pos1(i); p2 = D.pos2(i);
    r1 = get_ratio(c1,p1); r2 = get_ratio(c2,p2);
    line2(1e-6*[xon(c1)+p1 xon(c2)+p2],[yn(c1)+ys*r1 yn(c2)+ys*r2],'color',D.color{i},'linewidth',drlw);
  end
end

if ~isempty(P.inset_chr)
  % draw rearrangements affecting the gene:
  %   1. rearrangements within gene
  %   2. rearrangements from gene to other parts of its chromosome
  %   3. rearrangements from gene to the other chromosome

  D.inset1 = (D.chr1==igchr & D.pos1>=inset_start & D.pos1<=inset_end);
  D.inset2 = (D.chr2==igchr & D.pos2>=inset_start & D.pos2<=inset_end);
  for i=1:slength(D)
    if ~D.inset1(i) && ~D.inset2(i), continue; end
    c1 = find(D.chr1(i)==chrs); if isempty(c1), continue; end
    ratio1 = get_ratio(c1,D.pos1(i));
    if D.inset1(i)
      x1 = inset_left+isc*(D.pos1(i)-inset_start);
      if P.show_inset_copy_number
        y1 = gene_top+2*gene_height+ysi*ratio1;
      else
        y1 = gene_mid;
      end
    else
      x1 = 1e-6*(xon(c1)+D.pos1(i)); y1 = yn(c1)+ys*ratio1;
    end
    c2 = find(D.chr2(i)==chrs); if isempty(c2), continue; end
    ratio2 = get_ratio(c2,D.pos2(i));
    if D.inset2(i)
      x2 = inset_left+isc*(D.pos2(i)-inset_start);
      if P.show_inset_copy_number
        y2 = gene_top+2*gene_height+ysi*ratio2;
      else
        y2 = gene_mid;
      end
    else
      x2 = 1e-6*(xon(c2)+D.pos2(i)); y2 = yn(c2)+ys*ratio2;
    end
    line2([x1 x2],[y1 y2],'color',D.color{i},'linewidth',drlw);
  end

end

% set axes

xlim(1e-6*[min(leftn+as_column(xon)) max(rightn+as_column(xon))]);
ylim([ymin ymax]);

if P.no_labeling
  set(gca,'visible','off')
else
  fs = P.fontsize;
  xlabel('position (Mb)','fontsize',fs);

  params = {'fontsize',fs,'verticalalignment','top'};
  for c=1:nchrs
%    text(label_x,ytop(c),convert_chr_back(chrs(c)),params{:});
    text(label_x,yn(c) + 0.15 * (ytop(c)-yn(c)),convert_chr_back(chrs(c)),params{:});
  end

%  ylabel('copy ratio','fontsize',fs);

%  (DOESNT WORK ANYMORE:  due to multichr and/or inset?)
%  mn2 = md-0.9*compression_factor*(md-mn);
%  mx2 = md+0.9*compression_factor*(md-mn);
%  y = unique(round(((mn2:(mx2-mn2)/3:mx2))*10)/10);
%  yt = []; for c=nchrs:-1:1, yt = [yt yn(c)+ys*y]; end
%  ytl = text_to_lines(sprintf('%0.1f\n',repmat(y,1,nchrs)));
%  set(gca,'ytick',yt,'yticklabel',ytl,'fontsize',fs,'tickdir','out');

   set(gca,'ytick',[],'yticklabel',[]);
end
set(gcf,'color',[1 1 1]);



hold off

  function line2(x,y,varargin)
  % draw line indicating rearrangement
  %       (1) straight line if there is both X- and Y-offset
  %       (2) curved line if there is only negligible X-offset or only negligible Y-offset
  %       (3) dot if there is both negligible X- and Y-offset (extremely local)

    is_xoff = (abs(x(2)-x(1))>=1);
    is_yoff = (abs(y(2)-y(1))/y(2)>=0.1);
    if is_xoff && is_yoff
      line(x,y,varargin{:});
    elseif is_xoff || is_yoff 
      l2dir = -l2dir;
      xx = [x(1) (x(2)+x(1))/2 x(2)];
      yy = [y(1) y(1)*(1+l2dir) y(2)];
      cs = spline(xx,yy);
      xxx = linspace(x(1),x(2),10);
      plot(xxx,ppval(cs,xxx),varargin{:});
    else
      plot(x,y,'.-',varargin{:},'markersize',20);
    end
  end




  function r = get_ratio(c,pos)
    if strcmpi(P.line_anchors,'segseq_segs')
      idxx = find(segn{c}.left>pos,1);
      if isempty(idxx), idxx = slength(segn{c}); end
      r = segn{c}.ratios(max(1,idxx-1));
    elseif strcmpi(P.line_anchors,'segseq_dots')
       idxx = find(sn{c}.windows>pos,1);
      if isempty(idxx), idxx = length(sn{c}); end
      if idxx>1 && abs(sn{c}.ratios(idxx-1)-pos) < abs(sn{c}.ratios(idxx)-pos)
          idxx=idxx-1;
      end
      r = sn{c}.ratios(idxx);
    elseif strcmpi(P.line_anchors,'segfile')
      idxx = find(Sn{c}.start>pos,1);
      if isempty(idxx), idxx = slength(Sn{c}); end
      r = Sn{c}.ratio(max(1,idxx-1));
    end
  end

end    
