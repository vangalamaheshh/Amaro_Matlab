function hist_figure(fname,xx,th,hist_range)

if nargin==3
  hist_range=50;
end

lw=2;
figure(1); clf;
hist(xx(:),hist_range);
ax=axis;
if length(hist_range)~=1
  axis([ hist_range(1)-0.5*(hist_range(2)-hist_range(1)) ...
         hist_range(end)+0.5*(hist_range(2)-hist_range(1)) ax(3:4)]);
end
line([th th],[ax(3) ax(4)],'Color','r','Linewidth',lw);
print('-dpng','-r90','-f1',[fname '.png']);

