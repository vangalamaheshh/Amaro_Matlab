function fig = mutfig(inputfile,outputfile,P)
% mutfig(inputfile,outputfile)
%
% mutation figure-drawing script
%
% inputfile should have the following fields (columnns):
%
%   feature
%     protein (one per file)--gives characteristics of the overall construct
%     domain--one for each domain of the construct
%     mutation--one for each mutation
%     aabar (optional)--"black" or "grey"--prints "railroad-tracks"-style amino acid scalebar
%
%   type (for mutation: Missense, Nonsense, Indel, Splice)
%        (for domain: user-defined)
%        --> "type" determines color
%
%   start, end -- mutations are positioned at the average of these two
%
%   label -- (for mutations only)--text to display above protein
%
%   R,G,B -- color values between 0 and 1 (optional, defaults exist)
%
%   shape -- circle/square/triangle/inv_triangle (optional, default=circle)
%
% outputfile
%   filename of output file
%   must have extension of ".png", ".jpg", ".jpeg", ".tiff", ".tif", ".pdf", or ".eps"
% 
%%%   Amaro Taylor-Weiner w/ Mike Lawrence -- STAD TCGA AWG for Toshi

if ~exist('inputfile','var'), error('Please supply name of input file'); end
if ~exist(inputfile,'file'), error('Input file does not exist'); end

if ~exist('P','var') && exist('outputfile','var') && isstruct(outputfile)
  P=outputfile;
  clear outputfile
end

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'mutfig_style',1);

P.aa_scale = 'none';    % 'black'/'grey'/'none'
P.aa_text = 'none';
P.font_size = 17;

rand('twister',5489)     % so result is deterministic

%P = impose_default_value(P,'shape_meanings','hypermut_treat');
P = impose_default_value(P,'shape_meanings','none');
switch P.shape_meanings
  case 'hypermut_treat'
    shape_meanings = {'Non-treated';'Treated';'Hypermutated';''};
  case 'subtype'
    shape_meanings = {'Subtype 1';'Subtype 2';'Subtype 3';'Subtype 4'};
  case 'none'
    shape_meanings = {'','','',''};
  otherwise
    error('Unknown P.shape_meanings');
end

P = impose_default_value(P,'xpos','fit');   % 'fit' or 'actual'
P = impose_default_value(P,'labels','on');

% read input file
fprintf('Processing input file %s\n', inputfile);
params = []; params.lowercase_fieldnames = true;
F = load_struct(inputfile, params);
require_fields(F,{'feature';'type';'start';'end';'label'});
F.start = str2double(F.start);
F.end = str2double(F.end);
nf = length(F.feature);
F.count = repmat(1,nf,1);
if isfield(F,'direction') && isfield(F,'shape')
  idx = find(strcmp('down',F.direction) & strcmp('triangle',F.shape));
  F.shape(idx) = repmat({'inv_triangle'},length(idx),1);
end

% deal with colors
if all(isfield(F,{'r','g','b'}))
  F.r = str2double(F.r);
  F.g = str2double(F.g);
  F.b = str2double(F.b);
  F.color = cell(nf,1);
  for i=1:nf
    F.color{i} = [F.r(i) F.g(i) F.b(i)];
  end
else
  F.color = cell(nf,1);
end

% provide default/random colors where needed
for i=1:nf
  if isempty(F.color{i})
    if strcmpi('protein',F.feature{i})
      F.color{i} = [0.7 0.9 0.9];
    elseif strcmpi('mutation', F.feature{i})
      if strncmpi('missense', F.type{i},8)
        F.color{i} = [0 0.5 0.5];
      elseif strncmpi('indel', F.type{i},6)
        F.color{i} = [1 0.6 0];
      elseif strncmpi('splice', F.type{i},6)
        F.color{i} = [0.6 0 0.6];
      elseif strncmpi('nonsense', F.type{i},8)
        F.color{i} = [1 0 0];
      else % choose random bright color for this mutation type
        idx = grep(['^' F.type{i} '$'],F.type,1);
        tries=0;
        done=false;
        while(~done && tries<1000)
          tries=tries+1;
          c = 0.6*ones(1,3) + 0.4*rand(1,3);
          c(3*ceil(rand)) = 0;
          F.color{idx} = repmat(c,length(idx),1);
          done=true;
          for j=1:nf
            if ismember(j,idx), continue; end
            if isempty(F.color{j}), continue; end
            s=sum((F.color{i}-F.color{j}).^2);
            if s<0.08, done=false; end   % too similar to another color
          end
        end
      end
    elseif strcmpi('domain', F.feature{i})
      % choose selected colors for domains
      if strcmpi('Effector region',F.label{i})
          F.color{i}=[0 0 0];
      end
      if strcmpi('GTP',F.label{i})
          F.color{i}=[.3 .3 .3];
      end
      
      %idx = grep(['^' F.type{i} '$'],F.type,1);
      %tries=0;
      %done=false;
      %while(~done && tries<1000)
      %  tries=tries+1;
      %  c = 0.2*ones(1,3) + 0.6*rand(1,3);
      %  F.color{idx} = repmat(c,length(idx),1);
      %  done=true;
      %  for j=1:nf
      %    if ismember(j,idx), continue; end
      %    if isempty(F.color{j}), continue; end
      %    s=sum((F.color{i}-F.color{j}).^2);
      %    if s<0.08, done=false; end   % too similar to another color
       % end
      %end
    end
  end      
end

% default shape = circle
if ~isfield(F,'shape')
  F.shape = repmat({'circle'},nf,1);
end

% find protein characteristics
idx=find(strcmpi('protein',F.feature));
if length(idx)~=1, error('Definition must include exactly one "protein" record'); end
prot=[];
prot.color = F.color{idx};
prot.start = F.start(idx);
prot.end = F.end(idx);
prot.label = F.label{idx};
prot.length = prot.end-prot.start+1;

% see if there's an entry for AA scalebar
idx=find(strcmpi('aabar',F.feature));
for i=1:length(idx)
  if strcmpi(F.type{idx(i)},'black')
    P.aa_scale = 'black';
    P.aa_text = 'yes';
  elseif strcmpi(F.type{idx(i)},'grey')
    P.aa_scale = 'grey';
    P.aa_text = 'yes';
  end
end  

% define drawing area

% on width,height,left,bottom of screen display

if P.mutfig_style==1 || P.mutfig_style==2
  fig.width = 1200;
  fig.height = 600;
  fig.pos.left = 20;
  fig.pos.bottom = 100;
elseif P.mutfig_style==3
  fig.width = 1000;
  fig.height = 570;
  fig.pos.left = 20;
  fig.pos.bottom = 100;
end

fig.pos.top = fig.pos.bottom + fig.height;
fig.pos.right = fig.pos.left + fig.width;

p = get(gcf,'position');
pp = get(gcf, 'paperposition');
res = p(3)/pp(3);
set(gcf,'position',[fig.pos.left fig.pos.bottom fig.width fig.height]);
set(gcf,'paperposition',[0 0 fig.width fig.height]/res);
set(gcf,'papersize',[fig.width fig.height]/res);

set(gcf,'color',[1 1 1]);
fig.minx = 0;
fig.miny = 0;
fig.maxx = 1;
fig.maxy = 1;
set(gca,'position',[fig.minx fig.miny (fig.maxx-fig.minx) (fig.maxy-fig.miny)]);
set(gca,'xlim',[fig.minx fig.maxx]);
set(gca,'ylim',[fig.miny fig.maxy]);
axis off


if P.mutfig_style==1 || P.mutfig_style==2
  fig.prot.left = 0.1;
  fig.prot.width = 0.95-fig.prot.left;
  fig.prot.bottom = 0.375;
  fig.prot.height = 0.04;
elseif P.mutfig_style==3
  fig.prot.left = 0.26;
  fig.prot.width = 0.99-fig.prot.left;
  fig.prot.bottom = 0.25;
  fig.prot.height = 0.025;
end

fig.prot.right = fig.prot.left + fig.prot.width;
fig.prot.top = fig.prot.bottom + fig.prot.height;

fig.prot.xscale = fig.prot.width/prot.length;

fig.aatext.x = fig.prot.left * 0.8;
fig.aascale.bottom = fig.prot.bottom-0.025;
fig.aascale.height = 0.0133;
fig.aascale.tickheight = 0.0167;
fig.aascale.textspace = 0.0167;

% draw protein

rectangle('position',[fig.prot.left fig.prot.bottom fig.prot.width fig.prot.height],'FaceColor',prot.color);

% draw domains

D = reorder_struct(F,strcmpi('domain',F.feature));
for i=1:length(D.label)
  x1 = convx(D.start(i)-prot.start+1);
  x2 = convx(D.end(i)-prot.start+1);
  if D.end(i) >= prot.length, x2 = fig.prot.right; end
  w = x2-x1;
  y = fig.prot.bottom;
  h = fig.prot.height;
  rectangle('position',[x1 y w h], 'FaceColor', D.color{i});
end

% draw aa scale

if prot.length < 100, aax1=10; aax2=2;
elseif prot.length < 200, aax1=20; aax2=5;
elseif prot.length < 500, aax1=50; aax2=10;
elseif prot.length < 1000, aax1=100; aax2=20;
elseif prot.length < 2000, aax1=200; aax2=50;
elseif prot.length < 5000, aax1=500; aax2=100;
else aax1=1000; aax2=200; end

color = [0 0 0];
if strcmp(P.aa_scale,'grey'), color = [0.7 0.7 0.7]; end

if ~strcmp(P.aa_scale,'none')
  rectangle('position',[fig.prot.left fig.aascale.bottom fig.prot.width fig.aascale.height], ...
    'EdgeColor', color, 'FaceColor',color);
  for aa1 = prot.start:aax2*2:prot.end
    x1 = convx(aa1);
    x2 = convx(aa1+aax2);
    if (aa1+aax2) >= prot.length, x2 = fig.prot.right; end
    w = x2-x1;
    rectangle('position',[x1 fig.aascale.bottom w fig.aascale.height], ...
       'EdgeColor', color, 'FaceColor',[1 1 1]);
  end
  y = fig.aascale.bottom;
  fig.aatext.y = fig.aascale.bottom;
else
  y = fig.prot.bottom;
  fig.aatext.y = fig.prot.bottom - fig.aascale.tickheight;
end

if ~strcmp(P.aa_text,'none')
  text(fig.aatext.x,fig.aatext.y,'(AA)','fontsize',P.font_size,'color',[0 0 0],'rotation',90,...
     'verticalalignment','middle');
end
for aa = prot.start:aax1:prot.end
  x = convx(aa);
  line([x x], [y-fig.aascale.tickheight y],'color',color);
  text(convx(aa), y-fig.aascale.tickheight-fig.aascale.textspace, num2str(aa-1),...
      'fontsize',P.font_size,'color',color,'horizontalalignment','center');
end

% draw mutations
M = reorder_struct(F,strcmpi('mutation',F.feature));
if P.mutfig_style == 1
  % ORIGINAL GRAPHIC STYLE (June 2008)

  % get one record per mutation_name
  [tmp ui uj] = unique(M.label);
  M2 = reorder_struct(M,ui);
 
 
 
  for x=1:length(ui)
    M2.count(x) = sum(M.count(uj==x));
    M2.shape{x} = M.shape(uj==x);
    M2.color{x} = M.color(uj==x);
    M2.grouping{x}=M.grouping(uj==x);
  end
  [values,order] = sort([M2.start]);
 fields=fieldnames(M2);
  for i=1:length(fields)
     M2.(fields{i})=M2.(fields{i})(order);
 end
   
  most = max(M2.count);
  
  fig.symbol.xsize = 0.0117;
  fig.symbol.ysize = 0.0241;
  fig.symbol.yspacing = 0.01;
  
  fig.elbow.y1 = fig.prot.top + 0.0167;
  fig.elbow.y2 = fig.elbow.y1 + 0.0333;
  fig.elbow.y3 = fig.elbow.y2 + most*fig.symbol.ysize + (most+2)*fig.symbol.yspacing;
  fig.label.y = fig.elbow.y3 + 0.0167;
  
  % make more room in top of figure if necessary
  %if fig.height < fig.elbow.y3 + 0.2
    fig.maxy = fig.elbow.y3 + 0.2;
    set(gca,'YLim',[fig.miny fig.maxy]);
  %end

  M2.label = regexprep(M2.label,'_','\\_');    % protect underscores

  % increase size of asterisks
  M2.label = regexprep(M2.label,'\*',['\\fontsize\{' num2str(P.font_size*28/17) ...
                      '\}_\*\\fontsize\{' num2str(P.font_size) '\}']);
  
  % space mutation text
  nm = length(M2.label);
  a = (M2.end+M2.start)/2;
  [a ord] = sort(a);
  M2 = reorder_struct(M2,ord);
  th = text(a, ones(nm,1)*fig.label.y, M2.label,'rotation',90,'fontsize',P.font_size);
  s = prot.length / 85;
  
  switch P.xpos
   case 'fit'
    fprintf('Please ignore the following warning message:\n');
    x = find_text_pos(a,s,s,prot.start,prot.end);
   case 'actual'
    x = a;
   otherwise
    error('Unknown P.xpos');
  end
  
  for i=1:length(th)
    p = get(th(i),'position');
    p(1)=convx(x(i));
    set(th(i),'position',p);
    line(convx([a(i) a(i)]), [fig.prot.bottom fig.elbow.y1],'color',[0 0 0]);
    line(convx([a(i) x(i)]), [fig.elbow.y1 fig.elbow.y2],'color',[0 0 0]);
    if strcmp(P.labels,'on')
      line(convx([x(i) x(i)]), [fig.elbow.y2 fig.elbow.y3],'color',[0 0 0]);
  end
  
  y=fig.elbow.y2 + fig.symbol.yspacing/2 + fig.symbol.ysize/2;

  for k=1:M2.count(i);
    sx = convx(x(i))-fig.symbol.xsize/2;
    sy = y+fig.symbol.ysize/2;
    shape = M2.shape{i};
    if strcmp(M2.grouping{i}(k),'1')
    color = [67/255 170/255 224/255] ;
    end
    if  strcmp(M2.grouping{i}(k),'2')
        color = [238/255 71/255 35/255]
        
    end
    if  strcmp(M2.grouping{i}(k),'3')
    color = [146/255 128/255 186/255];
    end
    if  strcmp(M2.grouping{i}(k),'4')
    color = [125/255 194/255 66/255];
    end
   
    
    
    if strcmpi(shape,'circle')
      rectangle('position',[sx sy fig.symbol.xsize fig.symbol.ysize],...
                'curvature',[1 1],'facecolor',color,'edgecolor',[1 1 1]);
    elseif strcmpi(shape,'square')
      rectangle('position',[sx sy fig.symbol.xsize fig.symbol.ysize],...
                'curvature',[0 0],'facecolor',color,'edgecolor',[0 0 0]);
    elseif strcmpi(shape,'triangle')
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy sy sy+fig.symbol.ysize], color);
    elseif strcmpi(shape,'inv_triangle')
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy+fig.symbol.ysize sy+fig.symbol.ysize sy], color);
    else
      error('Unknown shape "%s"', shape);
    end
    y = y + (fig.symbol.ysize + fig.symbol.yspacing);
  end
  end
  
  if strcmp(P.labels,'off')
    delete(th);
  end

elseif P.mutfig_style==2
  % NEW STYLE (Jan 2013)
  gname = prot.label;
  text(fig.prot.left-0.005,fig.prot.bottom-0.005,gname,'fontsize',18-2*(length(gname)>8),...
       'fontsize',20,'horizontalalignment','right','verticalalignment','bottom');

  % mutation location window
  fig.mut = fig.prot;
  fig.mut.bottom = fig.prot.top + 0.01;
  fig.mut.height = 0.5;
  fig.mut.top = fig.mut.bottom + fig.mut.height;
  if ~isfield(M,'ypos'), M.ypos = rand(slength(M),1); end
  M = make_numeric(M,'ypos');

  fig.symbol.xsize = 0.0117;
  fig.symbol.ysize = 0.0241;

  % mutation locations
  for i=1:slength(M)
    x = convx((M.end(i)+M.start(i))/2);
    sx = x-fig.symbol.xsize/2;
    y = fig.mut.bottom + fig.mut.height * M.ypos(i);
    sy = y;
    shape = M.shape(i);
    
    color = M.color{i};
    if strcmpi(shape,'circle')
      rectangle('position',[sx sy fig.symbol.xsize fig.symbol.ysize],...
                'curvature',[1 1],'facecolor',color,'edgecolor',[0 0 0]);
    elseif strcmpi(shape,'square')
      rectangle('position',[sx sy fig.symbol.xsize fig.symbol.ysize],...
                'curvature',[0 0],'facecolor',color,'edgecolor',[0 0 0]);
    elseif strcmpi(shape,'triangle')
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy sy sy+fig.symbol.ysize], color);
    elseif strcmpi(shape,'inv_triangle')
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy+fig.symbol.ysize sy+fig.symbol.ysize sy], color);
    else
      error('Unknown shape "%s"', shape);
    end
  end

  fig.symbol.yspacing = 0.01;   % only needed for legend

elseif P.mutfig_style==3
  % NEWEST STYLE (May 2013): increase length:width ratio of main figure
  gname = prot.label;
  text(fig.prot.left-0.005,fig.prot.bottom-0.005,gname,'fontsize',18-2*(length(gname)>8),...
       'fontsize',20,'horizontalalignment','right','verticalalignment','bottom');

  % mutation location window
  fig.mut = fig.prot; fig.mut.height = 0.65;   % HEIGHT OF MUTATION DISPLAY AREA
  fig.mut.bottom = fig.prot.top + 0.01;
  fig.mut.top = fig.mut.bottom + fig.mut.height;
  if ~isfield(M,'ypos'), M.ypos = rand(slength(M),1); end
  M = make_numeric(M,'ypos');

  fig.symbol.xsize = 0.01;
  fig.symbol.ysize = 0.0175;

  % mutation locations
  for i=1:slength(M)
    x = convx((M.end(i)+M.start(i))/2);
    sx = x-fig.symbol.xsize/2;
    y = fig.mut.bottom + fig.mut.height * M.ypos(i);
    sy = y;
    shape = M.shape(i);
    color = M.color{i};
    if strcmpi(shape,'circle')
      rectangle('position',[sx sy 1.10*fig.symbol.xsize 1.10*fig.symbol.ysize],...
                'curvature',[1 1],'facecolor',color,'edgecolor',[0 0 0]);
    elseif strcmpi(shape,'square')
      rectangle('position',[sx sy fig.symbol.xsize fig.symbol.ysize],...
                'curvature',[0 0],'facecolor',color,'edgecolor',[0 0 0]);
    elseif strcmpi(shape,'triangle')
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy sy sy+fig.symbol.ysize], color);
    elseif strcmpi(shape,'inv_triangle')
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy+fig.symbol.ysize sy+fig.symbol.ysize sy], color);
    elseif strcmp(shape,'diamond')
      patch(sx+fig.symbol.xsize*1.20*[0 0.5 1 0.5], ...
            sy+fig.symbol.ysize*1.20*[0.5 0 0.5 1], color);
    else
      error('Unknown shape "%s"', shape);
    end
  end

  fig.symbol.yspacing = 0.01;   % only needed for legend


else
  error('unknown P.mutfig_style');
end

% draw legends

fig.legend.y = fig.aascale.bottom - fig.aascale.tickheight - fig.aascale.textspace - 0.086;
fig.legend.x1 = 0.1;
fig.legend.xspacing = 0.3;
fig.legend.symbol.xsize = fig.symbol.xsize;
fig.legend.symbol.ysize = fig.symbol.ysize;

if P.mutfig_style==1 || P.mutfig_style==2
  fig.legend.y = fig.aascale.bottom - fig.aascale.tickheight - fig.aascale.textspace - 0.086;
  fig.legend.symbol.yspacing = fig.symbol.yspacing * 1.5;
elseif P.mutfig_style==3
  fig.legend.y = fig.aascale.bottom - 0.05;
  fig.legend.symbol.yspacing = fig.symbol.yspacing * 0.90;
end

% domain legend

[tmp i j] = unique(D.type);
D2 = reorder_struct(D,sort(i));
D2.type = regexprep(D2.type,'_','\\_');   % protect underscores
x = fig.legend.x1;
y = fig.legend.y;
wmax = 0;
lenmax = 0;
for i=1:length(D2.type)
  rectangle('position', [x y fig.legend.symbol.xsize*4 fig.legend.symbol.ysize],...
     'edgecolor', [0 0 0], 'facecolor', D2.color{i});
  txt = D2.type{i};
  if length(txt)>50, txt=txt(1:50); end
  th = text(x+fig.legend.symbol.xsize*4+0.0167,y+fig.legend.symbol.ysize/2,txt,...
     'verticalalignment','middle','fontsize',P.font_size);
  e = get(th,'extent');
  wmax=max(wmax,e(3));
  lenmax=max(lenmax,length(txt));
  y = y - (fig.legend.symbol.ysize + fig.legend.symbol.yspacing);
end
ymin = y;

%x = x + fig.legend.xspacing;
%x = x + max(fig.legend.xspacing, wmax);   %<--- there is a problem with this
x = x + min(0.5,max([fig.legend.xspacing;wmax;0.012*lenmax]));

% shapes?

if ~ischar(shape_meanings) || ~strcmpi(shape_meanings,'none')
  shape_legend = false;
else
  shape_legend = true;
end

% MUTATION LEGEND

% choose most-saturated colors (where there is an ambiguity)
M.color_hsv = cell(slength(M),1);
M.color_sat = nan(slength(M),1);
for i=1:slength(M)
  M.color_hsv{i} = rgb2hsv(M.color{i});
  M.color_sat(i) = M.color_hsv{i}(2);
end
M = sort_struct(M,'color_sat',1);  % sort in forward order so that unique() picks the most-saturated (last) example

[tmp ui uj] = unique(M.groupchar);
M3 = reorder_struct(M,sort(ui,'descend'));
M3.groupchar = regexprep(M3.groupchar,'_','\\_');   % protect underscores
y = fig.legend.y;

    M3.color{4}=[67/255 170/255 224/255] ;
    M3.color{3}=[238/255 71/255 35/255];
    M3.color{2}=[146/255 128/255 186/255];
    M3.color{1}=[146/255 128/255 186/255];
M3.type=M3.groupchar;
for i=1:length(M3.groupchar)
  sx = x;
  sy = y;
  color = M3.color{i};
  if shape_legend
    shape = 'circle';
  else
    shape = M3.shape{i};
  end

  if strcmpi(shape,'circle')
    rectangle('position',[sx sy 1.10*fig.legend.symbol.xsize 1.10*fig.legend.symbol.ysize],...
              'curvature',[1 1],'facecolor',color,'edgecolor',[0 0 0]);
  elseif strcmpi(shape,'square')
    rectangle('position',[sx sy fig.legend.symbol.xsize fig.legend.symbol.ysize],...
              'curvature',[0 0],'facecolor',color,'edgecolor',[0 0 0]);
  elseif strcmpi(shape,'triangle')
    patch([sx sx+fig.legend.symbol.xsize sx+0.5*fig.legend.symbol.xsize], ...
          [sy sy sy+fig.legend.symbol.ysize], color);
  elseif strcmpi(shape,'inv_triangle')
    patch([sx sx+fig.legend.symbol.xsize sx+0.5*fig.legend.symbol.xsize], ...
          [sy+fig.legend.symbol.ysize sy+fig.legend.symbol.ysize sy], color);
  elseif strcmp(shape,'diamond')
    patch(sx+fig.legend.symbol.xsize*1.20*[0 0.5 1 0.5], ...
          sy+fig.legend.symbol.ysize*1.20*[0.5 0 0.5 1], color);
  else
    error('Unknown shape "%s"', shape);
  end

%  rectangle('position', [x y fig.legend.symbol.xsize fig.legend.symbol.ysize],...
%     'curvature', [1 1], 'edgecolor', [0 0 0], 'facecolor', M3.color{i});

  th = text(x+fig.legend.symbol.xsize+0.0167,y+fig.legend.symbol.ysize/2,M3.type{i},...
     'verticalalignment','middle','fontsize',P.font_size);
  e = get(th,'extent');
  wmax=max(wmax,e(3));
  y = y - (fig.legend.symbol.ysize + fig.legend.symbol.yspacing);
end
if y<ymin, ymin=y; end

x = x + fig.legend.xspacing;
%x = x + max(fig.legend.xspacing, wmax);   %<--- there is a problem with this

% symbol shape legend

if shape_legend
  y = fig.legend.y;
  grey = [0.7 0.7 0.7];
  white = [1 1 1];
  for i=1:4
    sy = y;
    sx = x;
    if isempty(shape_meanings{i}), continue; end
    if i==1 % circle
      rectangle('position',[sx sy 1.10*fig.symbol.xsize 1.10*fig.symbol.ysize],...
                'curvature',[1 1],'facecolor',white,'edgecolor',[0 0 0]);
    elseif i==2 % square
      rectangle('position',[sx sy fig.symbol.xsize fig.symbol.ysize],...
                'curvature',[0 0],'facecolor',white,'edgecolor',[0 0 0]);
    elseif i==3 % triangle
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy sy sy+fig.symbol.ysize], white);
    elseif i==4 % inv_triangle
      patch([sx sx+fig.symbol.xsize sx+0.5*fig.symbol.xsize], ...
            [sy+fig.symbol.ysize sy+fig.symbol.ysize sy], white);
    elseif i==5 % diamond
      patch(sx+fig.symbol.xsize*1.20*[0 0.5 1 0.5], ...
            sy+fig.symbol.ysize*1.20*[0.5 0 0.5 1], white);
    end
    th = text(x+fig.legend.symbol.xsize+0.0167,y+fig.legend.symbol.ysize/2,shape_meanings{i},...
              'verticalalignment','middle','fontsize',P.font_size);
    e = get(th,'extent');
    wmax=max(wmax,e(3));
    
    y = y - (fig.legend.symbol.ysize + fig.legend.symbol.yspacing);
  end
  if y<ymin, ymin=y; end
end


% make more room at bottom of figure if necessary
if ymin < 0.035
  fig.miny = ymin-0.035;
  fig.maxy = fig.maxy + 0.05;
  set(gcf, 'position',[fig.pos.left fig.pos.bottom fig.width fig.height + 0.15]);
  set(gca,'YLim',[fig.miny fig.maxy]);
end

if exist('outputfile','var') && ischar(outputfile) && ~isempty(outputfile)
  [P.output_format, P.output_res] = interpret_print_filename(outputfile);
  P.output_file = outputfile;
  fprintf('Generating output file %s\n', outputfile);
  device = ['-d' P.output_format];
  res = ['-r' num2str(P.output_res)];
  fname = P.output_file;
  print(device,res,fname);
  close all
end

fprintf('Done\n');


% helper function

  function x=convx(aan)
    x = fig.prot.left + (fig.prot.xscale.*(aan-1));
  end


end
