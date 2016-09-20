function slopes = qq(varargin)

idx = find(~cellfun(@isnumeric,varargin),1);
if isempty(idx)
  pv = varargin;
  vargs = {};
else
  pv = varargin(1:idx-1);
  vargs = varargin(idx:end);
end

if isodd(length(vargs))
  style = vargs{1};
  vargs = vargs(2:end);
else
  style = '.';
end

slopes = nan(length(pv),1);

tmp = {};
for i=1:length(pv)
  p = pv{i};
  if ~isnumeric(p)
    vargs = varargin(i:end);
    break
  end
  p = p(~isnan(p));
  p = max(0,p);
  q = calc_fdr_value(p);

  y=-log10(p);
  y(isinf(y)) = max(y(~isinf(y)));
  np = length(p);
  x=-log10((np:-1:1)'./(np+1));
  [y ord] = sort(y);
  p = p(ord);
  q = q(ord);
  % find slope at x=1.5
  pos=1.5;
  [zz xx] = min(abs(x-pos));
  if ~isempty(xx), slopes(i) = y(xx)/pos; end
  tmp = [tmp x y style];
  
  % find calc significance cutoff
  [zzz idx] = min(abs(q-0.1));
  threshx(i,1) = x(idx);
  threshy(i,1) = y(idx);
end

set(gcf,'DefaultAxesColorOrder',[...
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500    %very dark grey
    1         0.4       0         %orange
    0.6       0.6       1         %lt blue
    0.5       0         0.8       %violet
    0.35      0.35      0.35      %med dark gray
]);

plot(tmp{:},vargs{:});

xl=xlim; yl=ylim;
mn = min(xl(2),yl(2));
line([0 mn],[0 mn],'color',[0 0 0],'linestyle','--');
set(gca,'tickdir','out');

for i=1:length(pv)
  line(threshx(i)+0.07*[1 -1],threshy(i)*[1 1],'color',[1 0 0]);
end

if nargin==0
  clear slopes
end
