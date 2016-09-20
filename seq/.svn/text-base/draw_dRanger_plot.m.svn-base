function draw_dRanger_plot(X,P)
% draw_dRanger_plot(X,P)
%
% Mike Lawrence 2009-06-19

require_fields(X,{'chr1','pos1','str1','chr2','pos2','str2'});

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'TN_scale_factor',1);
P=impose_default_value(P,'dRanger_tumdir',[]);
P=impose_default_value(P,'dRanger_normdir',[]);
P=impose_default_value(P,'SegSeq_results_file',[]);
P=impose_default_value(P,'min_segment_depicted',10000);    % depicted segments are at least 1Mb
P=impose_default_value(P,'max_fractional_gap_depicted',0.2); % gaps larger than half the segment are suppressed
P=impose_default_value(P,'cnlog2ratio_mincap',-4);
P=impose_default_value(P,'cnlog2ratio_maxcap',4);
%P=impose_default_value(P,'tumreads_cap',30);
P=impose_default_value(P,'randseed',1234);

if isempty(P.dRanger_tumdir) | isempty(P.dRanger_normdir)
  error('Must specify dRanger_tumdir and dRanger_normdir');
end

if P.min_segment_depicted < 10000
  fprintf('Setting P.min_segment_depicted to 10,000 bp\n');
  P.min_segment_depicted = 10000;
end

try

% (1) build list of chromosomal regions
% (2) map them to the screen
% (3) draw them on the screen
% (4) draw the links

nx = slength(X);
all_chr = [X.chr1;X.chr2];
if any(isnan(all_chr) | all_chr<1 | all_chr>24) error('chr should be 1-24'); end
all_pos = [X.pos1;X.pos2];
if any(isnan(all_pos) | all_pos<1 | all_pos>300000000) error('pos should be 1-300M'); end

R = cell(24,1);   % Regions
numreg = zeros(24,1);
for c=1:24
  R{c}.start = [];
  R{c}.end = [];
  pos = all_pos(all_chr==c);
  while ~isempty(pos)
    mn = max(1,min(pos)-P.min_segment_depicted/2);
    mx = mn + P.min_segment_depicted - 1;
    while 1
%      fprintf('(%d) mn=%d mx=%d\n',c,mn,mx);
      incl = pos(pos<=mx);
%      fprintf('incl =');
%      for i=1:length(incl);fprintf(' %d',incl(i));end;fprintf('\n');
      if isempty(incl), break; end
      mx = max(incl) + P.min_segment_depicted/2 - 1;
      pos = setdiff(pos,incl);
    end
    R{c}.start = [R{c}.start; mn];
    R{c}.end = [R{c}.end; mx];
  end
  done = false;
  while ~done
    done = true;
    for i=1:slength(R{c})-1; j=i+1;
      totsize = R{c}.end(j) - R{c}.start(i);
      gapsize = R{c}.start(j) - R{c}.end(i);
      if gapsize/totsize < P.max_fractional_gap_depicted   % then merge
        done = false;
        R{c}.end(i) = R{c}.end(j);
        R{c}.end(j) = nan; R{c}.start(j) = nan;
      end
    end
    R{c}.end(isnan(R{c}.end))=[];
    R{c}.start(isnan(R{c}.start))=[];
  end
  R{c}.len = R{c}.end - R{c}.start + 1;
  numreg(c) = slength(R{c});
end
nchr = sum(numreg>0);

% assign screen real estate
% X = 0(left) to 1(right)
% Y = 0(bottom) to 1(top)

xleft = 0.1;
xright = 0.95;
ytop = 0.95;
ybottom = 0.05;

totalgapwidth = 0.1;
maxgapwidth = 0.1;

% for each chromosome determine:
% the total bp to be depicted
% the amount of horizontal real estate available
%   (after making room for gaps)
% the scale

bp_to_show = nan(24,1);
gapwidth = nan(24,1);
horiz_avail = nan(24,1);
for c=1:24
  bp_to_show(c) = sum(R{c}.len);
  if numreg(c)<2, gapwidth(c) = 0;
  else gapwidth(c) = totalgapwidth / (numreg(c)-1); end
  horiz_avail(c) = (xright-xleft)-gapwidth(c)*(numreg(c)-1);
end
scale = bp_to_show ./ horiz_avail;

% enforce uniform scale based on chromosome with the most to show
scale = max(scale);
horiz_used = bp_to_show / scale;

y1 = ytop - ((ytop-ybottom) * [0:nchr] / nchr);
ypos = (y1(1:end-1) + y1(2:end)) / 2;
row = 1;
for c=1:24, if numreg(c)==0, continue; end
  R{c}.ypos = ypos(row)*ones(numreg(c),1); row=row+1;
  if numreg(c)==1, gapsize = 0;
  else gapsize = min(maxgapwidth,((xright-xleft)-horiz_used(c))/(numreg(c)-1)); end
  tot_horiz = horiz_used(c) + gapsize * (numreg(c)-1);
  xpos = (xright+xleft)/2 - (tot_horiz/2);
  for r=1:numreg(c)
    R{c}.xstart(r,1) = xpos;
    R{c}.xend(r,1) = xpos + (R{c}.len(r) / scale);
    xpos = R{c}.xend(r,1) + gapsize;
  end
end

% load coverage

fprintf('Loading coverage data:');
C = cell(24,2);   % Coverage  T=1 N=2
for chr=1:24
  if numreg(chr)==0, continue; end
  fprintf(' chr%d',chr);
  for i=1:2
    if i==1, direc = P.dRanger_tumdir; else direc = P.dRanger_normdir; end
    f = [direc '/chr' num2str(chr) '.cov'];
    tmp = load_struct(f,'%f%f%f%f',0);
    C{chr,i} = rename_field(tmp,{'col1','col2','col3','col4'},...
      {'chr','start','end','count'});
end,end,fprintf('\n');

% extract coverage

fprintf('Extracting coverage for plot:\n');
for c=1:24
  R{c}.cnlog2ratio = cell(numreg(c),1);
  for r=1:numreg(c)
    fprintf('  chr%d region %d/%d: ',c,r,numreg(c));
    nkb = ceil(R{c}.len(r)/1000);
    R{c}.cnlog2ratio{r} = nan(nkb,1);
    for i=0:nkb-1
      if ~mod(i,1000), fprintf('%d/%d ',round(i/1000)+1,ceil(nkb/1000)); end
      pos = R{c}.start(r) + 1000*i;
      tidx = find(C{c,1}.start<=pos & C{c,1}.end>=pos);
      nidx = find(C{c,2}.start<=pos & C{c,2}.end>=pos);
      ttot = mean(C{c,1}.count(tidx));   %%%  NEEDS TO BE FIXED TO ACCOUNT FOR MISSING DATA
      ntot = mean(C{c,2}.count(nidx));
      cnratio = (ttot/P.TN_scale_factor) / ntot;
      cnlog2ratio = log2(cnratio);
      R{c}.cnlog2ratio{r}(i+1) = cnlog2ratio;
    end,fprintf('\n');
  end
end

% load SegSeq results if available

if ~isempty(P.SegSeq_results_file)
  fprintf('Loading SegSeq results.\n');
  SS = load_SegSeq_results(P.SegSeq_results_file);
  have_SS = true;
else
  have_SS = false;
end

%  draw figure
fprintf('Drawing figure...\n');
figure(1); clf; hold on
rectangle('position',[0 0 1 1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
set(gca,'position',[0.05 0.05 0.9 0.9],'xtick',[],'ytick',[]);

dna_slab = 0.04;

if have_SS
  y1 = 0;
  y2 = -dna_slab/2;
  yh = dna_slab/2;
else
  y1 = -dna_slab/2;
  yh = dna_slab;
end

for c=1:24, if numreg(c)==0, continue; end
  text(xleft/2,R{c}.ypos(1),sprintf('chr%d',c),'horizontalalign','center',...
    'verticalalign','middle');
  for r=1:numreg(c)
    st = R{c}.start(r); en = R{c}.end(r);
    xp = R{c}.xstart(r); yp = R{c}.ypos(r);
    % raw copy-number data
    n = length(R{c}.cnlog2ratio{r});
    for i=1:n
      cnlog2ratio = R{c}.cnlog2ratio{r}(i);
      if isnan(cnlog2ratio)
        k = [0.5 0.5 0.5];
      elseif cnlog2ratio>=0   % amp = red
        k = [1 1 1] - [0 1 1] * min(1,(cnlog2ratio/P.cnlog2ratio_maxcap));
      else    % del = blue
        k = [1 1 1] - [1 1 0] * min(1,(cnlog2ratio/P.cnlog2ratio_mincap));
      end
      rectangle('position',[xp+((i-1)*1000/scale) yp+y1 ...
        1000/scale yh],'facecolor',k,'edgecolor',k);
    end

    % SegSeq results
    if have_SS
      idx = find(SS.chr==c & SS.start<=en & SS.end>=st);
      for j=1:length(idx), s = idx(j);
        jst = max(st,SS.start(s));
        jen = min(en,SS.end(s));
        cnlog2ratio = log2(SS.copyratio(s));
        if isnan(cnlog2ratio)
          k = [0.5 0.5 0.5];
        elseif cnlog2ratio>=0   % amp = red
          k = [1 1 1] - [0 1 1] * min(1,(cnlog2ratio/P.cnlog2ratio_maxcap));
        else    % del = blue
          k = [1 1 1] - [1 1 0] * min(1,(cnlog2ratio/P.cnlog2ratio_mincap));
        end
%        fprintf('%d %d-%d %d [%.2f %.2f %.2f]\n',j,jst,jen,cnlog2ratio,k(1),k(2),k(3));
        rectangle('position',[xp+((jst-st)/scale) yp+y2 ...
          (jen-jst)/scale yh],'facecolor',k,'edgecolor',k);
      end
    end

    % black rectangle around whole segment
    line([xp xp R{c}.xend(r) R{c}.xend(r) xp],...
         yp+([-1 1 1 -1 -1]*(dna_slab/2)),'color',[0 0 0]);

    % draw position markers
    mark = ceil(st/1e6):floor(en/1e6); ndigits = 0;
    if isempty(mark), mark = (st+en)/2/1e6; ndigits = 2; end
    for i=1:length(mark)
      x = xp+((mark(i)*1e6-st)/scale);
      y = yp;
      line([x x],[y-dna_slab/2 y-1.2*dna_slab/2],'color',[0 0 0]);
      text(x,y-dna_slab,...
         sprintf(['%.' num2str(ndigits) 'f'],mark(i)),...
         'horizontalalign','center');
    end
%    text(xp,yp-dna_slab,sprintf('%d',round(st/1e6)),...
%      'horizontalalign','left');
%    text(R{c}.xend(r),yp-dna_slab,sprintf('%d',round(en/1e6)),...
%      'horizontalalign','right');
  end
end

% RANDOMIZE
rand('twister',P.randseed);

% draw rearrangements

curv=0.03;
for i=1:nx
  [x1 y1] = convert_coords(X.chr1(i),X.pos1(i));
  [x2 y2] = convert_coords(X.chr2(i),X.pos2(i));
  xcurv = min(curv, 0.5 * abs(x2-x1));
  ycurv = max(curv/3, 0.5 * abs(x2-x1));
  if X.str1(i)==0, xx1 = [x1 x1 x1+xcurv]; yy1 = [y1 y1 y1+ycurv];
              else xx1 = [x1 x1 x1-xcurv]; yy1 = [y1 y1 y1+ycurv]; end
  if X.str2(i)==0, xx2 = [x2+xcurv x2 x2]; yy2 = [y2+ycurv y2 y2];
              else xx2 = [x2-xcurv x2 x2]; yy2 = [y2+ycurv y2 y2]; end
  x = [x1 xx1 xx2 x2]; y = [y1 yy1 yy2 y2];
  c = spcrv([x;y]);
  linewidth = 2; % 5 * (min(P.tumreads_cap,X.tumreads(i))/P.tumreads_cap);
%  k = [0.8 0.8 0.8] - [0.8 0.8 0.8] * (min(P.tumreads_cap,X.tumreads(i))/P.tumreads_cap);
  k2 = rand(1,3);
%  plot(x1,y1,'o'); plot(x2,y2,'o');
%  plot(x,y,'o','color',k2);
  plot(c(1,:),c(2,:),'color',k2,'linewidth',linewidth);
end

xlim([0 1]);ylim([0 1]);
hold off

%keyboard

return

catch me; excuse(me); end

  function [x y] = convert_coords(chr,pos)
    x = nan;  y = nan;
    reg = find(pos>=R{chr}.start & pos<=R{chr}.end);
    if length(reg)==1
      y = R{chr}.ypos(reg);
      x = R{chr}.xstart(reg) + ((pos-R{chr}.start(reg))/scale);
    end
  end

end % main function
