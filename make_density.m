function [X,xr,yr]=make_density(ah,nx,ny)

if ~exist('ny','var')
  ny=nx;
end

axis(ah);
lh=findobj(ah,'Type','line','linestyle','none');
ax=axis;

xd=[];
yd=[];
for i=1:length(lh)
  xd=[xd get(lh(i),'XData')];
  yd=[yd get(lh(i),'YData')];
end

xr=ax(1):(ax(2)-ax(1))/(nx-1):ax(2);
yr=ax(3):(ax(4)-ax(3))/(ny-1):ax(4);
X=hist2d(xd,yd,xr,yr);
X=X(1:(end-1),1:(end-1));

