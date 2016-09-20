function draw_twochr_plot(S,seg,D,chra,chrb,P)
% OBSOLETE: replaced by draw_multichr_plot.m

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'draw_grey_rects',false);
P = impose_default_value(P,'segfile',[]);
P = impose_default_value(P,'inset_gene',[]);
P = impose_default_value(P,'refseq_entry',[]);
P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'midstrip',0.1+0.3*~isempty(P.inset_gene));
P = impose_default_value(P,'no_labeling',false);
P = impose_default_value(P,'show_segseq_dots',true);
P = impose_default_value(P,'show_segseq_segs',false);
P = impose_default_value(P,'line_anchors','segfile');
P = impose_default_value(P,'expanded_color_scheme',false);
P = impose_default_value(P,'xoa',0);
P = impose_default_value(P,'xob',0);
P = impose_default_value(P,'compression_factor',[]);
P = impose_default_value(P,'fontsize',16);

sa = reorder_struct(S.RATIOS,S.RATIOS.chr==chra);
sb = reorder_struct(S.RATIOS,S.RATIOS.chr==chrb);
sega = reorder_struct(seg.SEG,seg.SEG.chr==chra);
segb = reorder_struct(seg.SEG,seg.SEG.chr==chrb);

if ~isempty(P.inset_gene)
  if isempty(P.refseq_entry)
    R = load_refseq(P.build);
    Rc = reorder_struct(R,strcmpi(R.gene,P.inset_gene));
    if slength(Rc)==0, error('%s not found in RefSeq',P.inset_gene); end
    [tmp ord] = sort(Rc.n_exons,'descend');
    Rc = reorder_struct(Rc,ord(1));
  else
    Rc = P.refseq_entry;
  end
  inset_gene_chr = convert_chr(Rc.chr);
  if inset_gene_chr~=chra && inset_gene_chr~=chrb, error('%s is not on chr%d or chr%d',P.inset_gene,chra,chrb); end
end

% copy-number dynamic range
v = [sa.ratios;sb.ratios];
mn = min(v);
mx = max(v);
md = (mn+mx)/2;

%  chromosomal coordinates range
pa = double([sa.windows;sega.left;sega.right]);
pb = double([sb.windows;segb.left;segb.right]);
lefta = min(pa);
leftb = min(pb);
righta = max(pa);
rightb = max(pb);
left = min(lefta,leftb);
right = max(righta,rightb);
width = right-left+1;

% SET UP DRAWING AREAS

ymin = 0; ymax = 1;
topstrip = 0.05;
btmstrip = 0.05;
midstrip = P.midstrip; %-0.1; %0.1;

chrstrip = (1-(topstrip+btmstrip+midstrip))/2;
if (chrstrip<=0) error('margin strips too big'); end

if (mx-mn)>=2, compression_factor=1.2; else compression_factor = 1.2/(mx-mn); end

if ~isempty(P.compression_factor), compression_factor = P.compression_factor; end

ys = (0.5*chrstrip) / (compression_factor*(md-mn));

ya = 1 - topstrip - 0.5*chrstrip - (md*ys);
yb = 0 + btmstrip + 0.5*chrstrip - (md*ys);

% x-offsets
xoa = P.xoa;
xob = P.xob;
label_x = double(width*0.01*1e-6);

% DRAW PLOT

clf;hold on

if P.draw_grey_rects
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

% draw segfile segments (if available)
if ~isempty(P.segfile)
  S = load_cnseg_file(P.segfile);
  S.copyratio = 2.^S.segmean;
  min_nprobes = 8;
  S = reorder_struct(S,S.nprobes>=min_nprobes);
  Sa = reorder_struct(S,S.chr==chra);
  for i=1:slength(Sa)
    line(1e-6*(xoa+[Sa.start(i) Sa.end(i)]),ya+ys*Sa.copyratio(i)*[1 1],'color',[1 0.6 0],'linewidth',4);
  end
  Sb = reorder_struct(S,S.chr==chrb);
  for i=1:slength(Sb)
    line(1e-6*(xob+[Sb.start(i) Sb.end(i)]),yb+ys*Sb.copyratio(i)*[1 1],'color',[1 0.6 0],'linewidth',4);
  end
end


% segseq
if P.show_segseq_dots
  scatter(1e-6*(xoa+sa.windows),ya+ys*sa.ratios,'.k')
  scatter(1e-6*(xob+sb.windows),yb+ys*sb.ratios,'.k')
end

if P.show_segseq_segs
  for i=1:slength(sega)
    line(1e-6*(xoa+[sega.left(i),sega.right(i)]),[1 1]*(ya+ys*sega.ratios(i)),'color','green','linewidth',2.5);
  end
  for i=1:slength(segb)
    line(1e-6*(xob+[segb.left(i),segb.right(i)]),[1 1]*(yb+ys*segb.ratios(i)),'color','green','linewidth',2.5); 
  end
end

drlw = 1.5;   % dRanger linewidth
intercolor = [0.8 0 0.8];
intracolor = [0.2 0.9 0.2];

% dRanger intra class colors
classes = {'deletion','inversion','long_range','tandem_dup'};
if P.expanded_color_scheme
  class_colors = {[0 0 1],[0 1 1],[0.2 0.9 0.2],[1 0 0]};
else
  class_colors = repmat({intracolor},1,4);
end
D.color = map_across(D.class,classes,class_colors);

if isempty(P.inset_gene)
  % if no inset gene, then show all rearrangements between the two chromosomes
 
  % dRanger inter
  idx = find(D.chr1==chra & D.chr2==chrb);
  for i=1:length(idx)
    p1 = D.pos1(idx(i)); p2 = D.pos2(idx(i));   
    r1 = get_ratio('a',p1); r2 = get_ratio('b',p2);
    line(1e-6*[xoa+p1 xob+p2],[ya+ys*r1 yb+ys*r2],'color',intercolor,'linewidth',drlw);
  end
  idx = find(D.chr1==chrb & D.chr2==chra);
  for i=1:length(idx)
    p1 = D.pos1(idx(i)); p2 = D.pos2(idx(i));
    r1 = get_ratio('a',p1); r2 = get_ratio('b',p2);
    line(1e-6*[xob+p1 xoa+p2],[yb+ys*r1 ya+ys*r2],'color',intercolor,'linewidth',drlw);
  end

  % dRanger intra  
  idx = find(D.chr1==chra & D.chr2==chra);
  for i=1:length(idx)
    p1 = D.pos1(idx(i)); p2 = D.pos2(idx(i));
    r1 = get_ratio('a',p1); r2 = get_ratio('a',p2);
    line(1e-6*[xoa+p1 xoa+p2],[ya+ys*r1 ya+ys*r2],'color',D.color{idx(i)},'linewidth',drlw);
  end
  idx = find(D.chr1==chrb & D.chr2==chrb);
  for i=1:length(idx)
    p1 = D.pos1(idx(i)); p2 = D.pos2(idx(i));
    r1 = get_ratio('b',p1); r2 = get_ratio('b',p2);
    line(1e-6*[xob+p1 xob+p2],[yb+ys*r1 yb+ys*r2],'color',intracolor,'linewidth',drlw);
  end

else   % ~isempty(P.inset_gene)

  % set up gene drawing area
  gene_horiz_margin = 0.15;
  gene_left = 1e-6*(left+width*gene_horiz_margin); 
  gene_right = 1e-6*(right-width*gene_horiz_margin);
  gene_width = gene_right - gene_left;

  gene_start = Rc.tx_start;
  gene_end = Rc.tx_end;
  gene_len = gene_end - gene_start + 1;
  gsc = gene_width / gene_len;

  gene_height = 0.08;
  gene_mid = btmstrip+chrstrip+(0.5*midstrip);
  gene_top = gene_mid+(0.5*gene_height);
  gene_btm = gene_mid-(0.5*gene_height); 

  % draw intron/exon diagram
  rectangle('position',[gene_left gene_btm gene_width gene_height],'facecolor',[0.8 0.8 0.8]);
  text(label_x,gene_mid,P.inset_gene,'verticalalignment','middle','fontsize',P.fontsize);  % show gene name
  for i=1:Rc.n_exons
    st = Rc.exon_starts{1}(i); en = Rc.exon_ends{1}(i);
    rectangle('position',[gene_left+gsc*(st-gene_start) gene_btm gsc*(en-st) gene_height],'facecolor',[0 0 0]);
  end

  % draw rearrangements affecting the gene:
  %   1. rearrangements within gene
  %   2. rearrangements from gene to other parts of its chromosome
  %   3. rearrangements from gene to the other chromosome

  D.inset1 = (D.chr1==inset_gene_chr & D.pos1>=gene_start & D.pos1<=gene_end);
  D.inset2 = (D.chr2==inset_gene_chr & D.pos2>=gene_start & D.pos2<=gene_end);

  % rearrangements from gene to other parts of its chromosome
  for i=1:slength(D)
    if D.inset1 & D.inset2    % rearrangement within gene
      % (none yet)
    elseif D.inset1
      if D.chr2 == inset_gene_chr   % rearrangement from gene to another part of its chromosome

      else   % rearrangement from gene to the other chromosome
      end 
    elseif D.inset2
      if D.chr2 == inset_gene_chr   % rearrangement from gene to another part of its chromosome

      else   % rearrangement from gene to the other chromosome
      end
    end
  end

%    p1 = D.pos1(idx(i)); x1 = gene_left+gsc*(p1-gene_start);
%    p2 = D.pos2(idx(i));
%    r1 = get_ratio('b',p1); r2 = get_ratio('b',p2);
%    line(1e-6*[xob+p1 xob+p2],[yb+ys*r1 yb+ys*r2],'color',intracolor,'linewidth',drlw);


%  for i=1:slength(Dc)
%    if Dc.chr2(i)==5, col = [0.2 0.8 0.2]; else col = [0.8 0.2 0.8]; end
%    line([1 1]*Dc.pos1(i),[0.7 1],'color',col,'linewidth',2.5);
%  end

end

% set axes

xlim(1e-6*[left right]);
ylim([ymin ymax]);

if P.no_labeling
  set(gca,'visible','off')
else
  fs = P.fontsize;
  xlabel('position (Mb)','fontsize',fs);

  params = {'fontsize',fs,'verticalalignment','top'};
  text(label_x,btmstrip+midstrip+2*chrstrip+0.02,convert_chr_back(chra),params{:});
  text(label_x,btmstrip+chrstrip+0.02,convert_chr_back(chrb),params{:});
  
  ylabel('copy ratio','fontsize',fs);
  mn2 = md-0.9*compression_factor*(md-mn);
  mx2 = md+0.9*compression_factor*(md-mn);
  y = unique(round(((mn2:(mx2-mn2)/3:mx2))*10)/10);
  yt = [yb+ys*y ya+ys*y];
  ytl = text_to_lines(sprintf('%0.1f\n',[y y]));
  set(gca,'ytick',yt,'yticklabel',ytl,'fontsize',fs,'tickdir','out');
end
set(gcf,'color',[1 1 1]);



hold off


  function r = get_seg_ratio_at_pos(seg,pos)
    idxx = find(seg.left<=pos & seg.right>=pos);
    if isempty(idxx), r = nan; else r = seg.ratios(idxx(1)); end
  end

  function r = get_ratio_at_pos(ratio,pos)
    idxx = find(ratio.windows>pos,1);
    if isempty(idxx), idxx = length(ratio); end
    r = ratio.ratios(idxx);
  end

  function r = get_ratio_in_segfile(segfile,pos)
    idxx = find(segfile.start<=pos & segfile.end>=pos);
    if isempty(idxx), r = nan; else r = segfile.ratio(idxx); end
  end

  function r = get_ratio(ab,pos)
    if strcmpi(P.line_anchors,'segseq_segs')
      if strcmp(ab,'a'), r = get_seg_ratio_at_pos(sega,pos);
      elseif strcmp(ab,'b'), r = get_seg_ratio_at_pos(segb,pos); end
    elseif strcmpi(P.line_anchors,'segseq_dots')
      if strcmp(ab,'a'), r = get_ratio_at_pos(sa,pos);
      elseif strcmp(ab,'b'), r = get_ratio_at_pos(sb,pos); end
%      r1 = get_ratio_at_pos(sb,p1);
%      r2 = get_ratio_at_pos(sb,p2);
    elseif strcmpi(P.line_anchors,'segfile')
      if strcmp(ab,'a'), r = get_ratio_in_segfile(Sa,pos);
      elseif strcmp(ab,'b'), r = get_ratio_in_segfile(Sb,pos); end
%      r1 = get_ratio_in_segfile(Sb,p1);
%      r2 = get_ratio_in_segfile(Sb,p2);
    end
  end

end    
