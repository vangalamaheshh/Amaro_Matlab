function pout = ebadjust(pin)

y=-log10(pin);
y(isinf(y)) = max(y(~isinf(y)));
n = length(pin);
x=-log10((1:n)'./(n+1));
[x ordx] = sort(x);
[y ordy] = sort(y);

% FIT
% iterate over possible windows
%    bounded on left by end of run of 1's
%    bounded on right by where y passes 1/n


% (exclude initial run of 1's)

leftbound = x(find(y>0,1));
rightbound = x(find(y>-log10(1/n)),1);
if isempty(rightbound)
  rightbound = -log10(1/n);
end

step=0.05;
best=[]; best.score = -inf;
for left=leftbound:step:rightbound-step
  for right=left+step:step:rightbound
    idx = find(x>=left & x<=right);
    if isempty(idx), continue; end
    x1=x(idx(1)); x2=x(idx(end));
    y1=y(idx(1)); y2=y(idx(end));
    if (x2-x1)<0.3, continue; end  % require fits to span at least 1/2 order of magnitude
    a = (x1-x2)/(y1-y2);
    if (a>3 || a<(1/3)), continue; end % don't allow extreme slopes
    b = x1-a*y1;
    if (abs(b)>5), continue; end  % don't allow extreme offsets
    yfit = a*y(idx)+b;
    resid = abs(x(idx)-yfit);
    dx = diff(x(idx));
    auc = sum(resid(1:end-1).*dx);
%    keyboard
%    closeness = (1./resid);
%    closeness(resid>0.5) = 0;
%    closeness(resid<0.1) = 10;
%    score = mean(closeness)*(x2-x1);
    score = (x2-x1).^2/auc;
    if score>best.score
      best.a = a;
      best.b = b;
      best.x1 = x1;
      best.x2 = x2;
      best.y1 = y1;
      best.y2 = y2;
      best.left = left;
      best.right = right;
      best.score = score;
end,end,end

clf,hold on
plot(x,y,'b.');

if ~isfield(best,'a')
  fprintf('WARNING: failed to fit\n');
  yfit = y;
else
  yfit = (best.a*y)+best.b;
  plot(x,yfit,'r.');
  plot(best.x1,best.x1,'r.','markersize',35);
  plot(best.x2,best.x2,'r.','markersize',35);
  plot(best.x1,best.y1,'b.','markersize',35);
  plot(best.x2,best.y2,'b.','markersize',35);
  xl=xlim; yl=ylim;
  x1=best.x1-0.5; x2=best.x2+0.5;
  slope = (best.y2-best.y1)/(best.x2-best.x1);
  y1=best.y1-0.5*slope; y2=best.y2+0.5*slope;
  line([x1 x2],[y1 y2],'color',[0 0 0],'linestyle','-');
  xlim(xl); ylim(yl);
end

xl=xlim; yl=ylim;
mn = min(xl(2),yl(2));
line([0 mn],[0 mn],'color',[0 0 0],'linestyle','--');
hold off
%png

% apply transformation to p-values
yout = nan(n,1);
yout(ordy) = yfit;
pout = 10.^(-yout);
pout(pout>1) = 1;
pout(pin==1) = 1;
pout(pin==0) = 0;
