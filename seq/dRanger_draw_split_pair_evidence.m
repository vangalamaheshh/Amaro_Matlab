function dRanger_draw_split_pair_evidence(R,W,widx,ridx,P,Rf)
% dRanger_draw_split_pair_evidence(R,sidx,ridx,P)
%
%   R = dRanger results struct
%   W = weird pairs matrix (id chr1 str1 st1 en1 chr2 str2 st2 en2 good1 good2 readgroup switch)
%   widx = mapping from W to R
%   ridx = which entry of R to illustrate
%   P = parameters
%   Rf = Refseq from load_refseq
%
% Mike Lawrence 2010-11-27

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'randseed',54);
P = impose_default_value(P,'max_pairs_to_show',1000);
P = impose_default_value(P,'read_min_display_width_pct',0.5);
P = impose_default_value(P,'draw_genes',true);
P = impose_default_value(P,'genes_to_suppress',{});
P = impose_default_value(P,'fontsize',15);
P = impose_default_value(P,'gene_names_fontsize',P.fontsize);
P = impose_default_value(P,'format_for','screen');
if strcmpi(P.format_for,'screen')
  P = impose_default_value(P,'read_linewidth',6);
  P = impose_default_value(P,'fragment_linewidth',2);
  P = impose_default_value(P,'gca_visible','on');
  P = impose_default_value(P,'gene_box_color',[0.8 0.8 0.8]);
  P = impose_default_value(P,'gene_strip_gapsize',0.03);
elseif strcmpi(P.format_for,'pdf')
  P = impose_default_value(P,'read_linewidth',3);
  P = impose_default_value(P,'fragment_linewidth',1);
  P = impose_default_value(P,'gca_visible','off');
  P = impose_default_value(P,'gene_box_color',[0.7 0.7 0.75]);
  P = impose_default_value(P,'gene_strip_gapsize',0.04);
else
  error('unknown P.format_for');
end

demand_fields(R,{'num','chr1','chr2','str1','str2','min1','min2','max1','max2'});
if size(W,2)~=13, error('W should have 13 columns'); end

% rearrangement info
Side{1}.chr = R.chr1(ridx);
Side{2}.chr = R.chr2(ridx);
Side{1}.str = R.str1(ridx);
Side{2}.str = R.str2(ridx);
Side{1}.min = R.min1(ridx);
Side{2}.min = R.min2(ridx);
Side{1}.max = R.max1(ridx);
Side{2}.max = R.max2(ridx);

% extract pairs
w = W(widx==R.num(ridx),2:9);
nw = size(w,1);

figure(1)
scatter(w(:,3),w(:,7))

% reverse end1/end2 if more convenient
if Side{1}.str==1 && Side{2}.str==0
  tmp=Side{1}; Side{1}=Side{2}; Side{2}=tmp;
  w = w(:,[5:8 1:4]);
end

Side{1}.reverse = (Side{1}.str==1);
Side{2}.reverse = (Side{2}.str==0);

% order pairs randomly
rand('twister',P.randseed);
w = w(randperm(nw),:);

if ~isempty(P.max_pairs_to_show) && nw>P.max_pairs_to_show
  nw = P.max_pairs_to_show;
  w = w(1:nw,:);
end

% choose genome ranges to show; compute screen layout
marg1 = 0;%round(0.1*(Side{1}.max-Side{1}.min));
marg2 = 0;%round(0.1*(Side{2}.max-Side{2}.min));
Side{1}.left = Side{1}.min-marg1; Side{1}.right = Side{1}.max+marg1;
Side{2}.left = Side{2}.min-marg2; Side{2}.right = Side{2}.max+marg2;

%leftgap = marg1*2;
%rightgap = marg2*2;
%midgap = leftgap+rightgap;

m1 = round(0.1*(Side{1}.max-Side{1}.min));
m2 = round(0.1*(Side{2}.max-Side{2}.min));
leftgap = m1*2;
rightgap = m2*2;
midgap = leftgap+rightgap;

Side{1}.start = leftgap; Side{1}.width = Side{1}.right-Side{1}.left;
Side{2}.start = leftgap + Side{1}.width + midgap; Side{2}.width = Side{2}.right - Side{2}.left;
totwidth = Side{1}.width+Side{2}.width+leftgap+midgap+rightgap;
Side{1}.end = Side{1}.start+Side{1}.width;
Side{2}.end = Side{2}.start+Side{2}.width;

reads_strip = nw*2+2;
if P.draw_genes
  gene_strip_top_gap = P.gene_strip_gapsize*reads_strip;
  gene_strip = 0.03 * reads_strip;
  gene_strip_btm_gap = P.gene_strip_gapsize*reads_strip;
  ymin = -(gene_strip_top_gap+gene_strip+gene_strip_btm_gap);
else
  ymin = 0;
end
ymax = reads_strip;

% draw figure
figure(2)
clf;hold on;
minreadwidth = totwidth*P.read_min_display_width_pct/100;
grey = ones(1,3)*0.8;
rectangle('position',[Side{1}.end ymin-1 midgap (ymax-ymin)+10],'facecolor',grey);

% draw reads
for i=1:nw
  % read1
  st = w(i,3); en = w(i,4);
  if en-st<minreadwidth
    mid = (en+st)/2;
    st = mid-(minreadwidth/2);    
    en = mid+(minreadwidth/2);
  end
  if ~Side{1}.reverse
    lf = Side{1}.start+(st-Side{1}.left);
    rt = Side{1}.start+(en-Side{1}.left);
  else
    lf = Side{1}.end-(en-Side{1}.left);
    rt = Side{1}.end-(st-Side{1}.left);
  end
  line([lf rt],i*2*[1 1],'linewidth',P.read_linewidth,'color',[0.2 0.2 0.8]);
  line([rt Side{1}.end],i*2*[1 1],'linewidth',P.fragment_linewidth,'color',[0 0 0]);

  % line connecting the two reads
  line([Side{1}.end Side{2}.start],i*2*[1 1],'linewidth',P.fragment_linewidth,'color',[0 0 0],'linestyle','-');

  % read2
  st = w(i,7); en = w(i,8);
  if en-st<minreadwidth
    mid = (en+st)/2;
    st = mid-(minreadwidth/2);
    en = mid+(minreadwidth/2);
  end
  if ~Side{2}.reverse
    lf = Side{2}.start+(st-Side{2}.left);
    rt = Side{2}.start+(en-Side{2}.left);
  else
    lf = Side{2}.end-(en-Side{2}.left);
    rt = Side{2}.end-(st-Side{2}.left);
  end
  line([Side{2}.start lf],i*2*[1 1],'linewidth',P.fragment_linewidth,'color',[0 0 0]);
  line([lf rt],i*2*[1 1],'linewidth',P.read_linewidth,'color',[0.8 0.2 0.2]);
end % next read

if ~Side{1}.reverse
  xtl = [Side{1}.left Side{1}.right];
else
  xtl = [Side{1}.right Side{1}.left];
end
if ~Side{2}.reverse
  xtl = [xtl Side{2}.left Side{2}.right];
else
  xtl = [xtl Side{2}.right Side{2}.left];
end  
xt = [Side{1}.start Side{1}.start+Side{1}.width Side{2}.start Side{2}.start+Side{2}.width];
set(gca,'xtick',xt,'xticklabel',num2cellstr(xtl),'tickdir','out','ytick',[]);

% draw genes (if specified)
if P.draw_genes
  if ~exist('Rf','var')
    Rf = load_refseq(P.build); Rf.chr = convert_chr(Rf.chr);
  end

  for side=1:2
    G = reorder_struct(Rf,Rf.chr==Side{side}.chr & Rf.tx_start<=Side{side}.right & Rf.tx_end>=Side{side}.left);
    G = reorder_struct(G,~ismember(G.gene,P.genes_to_suppress));
    [u ui uj] = unique(G.gene);
    G = reorder_struct(G,ui);
    G = sort_struct(G,'tx_start');
    toggle = false;
    for gi=1:slength(G)
      % gene box
      [x y] = calculate_jag_box(G.tx_start(gi),G.tx_end(gi),Side{side}.left,Side{side}.right);
      if ~Side{side}.reverse
        px = Side{side}.start+(x-Side{side}.left);
      else
        px = Side{side}.end-(x-Side{side}.left);
      end
      py = -(gene_strip_top_gap+gene_strip)+y*gene_strip;
      patch(px,py,P.gene_box_color);
      % gene name
      tx = (px(end)+px(end-1))/2;
      plusstrand = strcmp(G.strand{gi},'+');
      glabel = make_arrow_label(G.gene{gi},xor(~plusstrand,Side{side}.reverse));
      if toggle
        ty = 0;
      else
        ty = -(gene_strip_top_gap+gene_strip);
      end
      ty = ty-(gene_strip*0.15);
      toggle = ~toggle;
      params = {'verticalalignment','top','fontsize',P.gene_names_fontsize,'horizontalalignment','center'};
      text(tx,ty,glabel,params{:});
      % gene exons
      for ei=1:G.n_exons(gi)
        st = G.exon_starts{gi}(ei); en = G.exon_ends{gi}(ei);
        if (st<Side{side}.right && en>Side{side}.left)
          [x y] = calculate_jag_box(st,en,Side{side}.left,Side{side}.right);
          if ~Side{side}.reverse
            px = Side{side}.start+(x-Side{side}.left);
          else
            px = Side{side}.end-(x-Side{side}.left);
          end
          py = -(gene_strip_top_gap+gene_strip)+y*gene_strip;
          patch(px,py,[0 0 0]);
        end
      end % next exon
    end % next gene
  end % next side 
end


% chromosome labels
ty = 0-ymax*0.12;
params = {'fontsize',P.fontsize,'verticalalignment','top','horizontalalignment','center'};
for side=1:2
  chrlabel = make_arrow_label(decell(convert_chr_back(Side{side}.chr)),Side{side}.reverse);
  text(Side{side}.start+Side{side}.width/2,ty,chrlabel,params{:});
end



ylim([ymin ymax]);
%line(Side{1}.end+midgap/2*[1 1],[ty ymax],'


set(gca,'visible',P.gca_visible);
% set(gcf,'color',[1 1 1]);
hold off

set(gcf,'position',[83 103 1150 839]);


  function [x y] = calculate_jag_box(st,en,winleft,winright)
    jag = 0.005*totwidth;   % width of 'jagged cut' for truncated boxes
    if st>=winleft
      x = [st st];
      y = [0 1];
    else
      x = winleft-jag*[0 1 0 1 0 1 0 1 0];
      y = [0 0:6 6]/6;
    end
    if en<=winright
      x = [x en en];
      y = [y 1 0];
    else
      x = [x winright+jag*[0 1 0 1 0 1 0]];
      y = [y [6:-1:0]/6];
    end
    x = [x x(1)];
    y = [y y(1)];
  end

  function label = make_arrow_label(base_label, reverse_flag)
    if ~reverse_flag
%      label = ['--' base_label '-->'];
      label = [base_label '\rightarrow'];
    else
%      label = ['<--' base_label '--'];
      label = ['\leftarrow' base_label];
    end
  end

end  % main function
