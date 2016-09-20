function [h,xb]=myhistT(x,b,t)
%function [h,xbin,Stats]=myhistT(x,b,t,arg)
% plots two histograms for x in bins b as stairstepped lines 
% for condition t (true or false)
%
% inputs: x: vector of values to histgram
%         b: bins (same convention as hist)
%         t: vector of condition of same length as x
%
% outputs:  h: handle to stairs
%           xbin: bin centers vector
%
if (isnumeric(t))
    t=t~=0;
end
[n1,xb]=hist(x(t),b); 
[n2,xb]=hist(x(~t),xb); 
[q1,n1]=hist2stair(xb,n1);
[q2,n2]=hist2stair(xb,n2);
h=stairs(q1,[n1 n2 n2*0]); 
set(h,'linewidth',2);  
set(h(3),'color','k')
%set(gca,'xscale','log')
%    set(gca,'xlim',minmax(x),'xtick',10.^(1:6))
%    text('units','normalized','position',[0.05 0.95],'string',cellstr(A(a1)),'fontsize',14, ...
%    'HorizontalAlignment','left','VerticalAlignment','top','interpreter','tex','backgroundcolor','w');

function [x,y]=hist2stair(x1,y1)

x1=x1(:)';
N=size(y1);
if (N(1)==length(x1)) 
   y1=y1';
end
% bin widths
dx=diff(x1); 
% tack on a dx at the end
dx = [dx dx(end)];
% move x to bin starts
x=x1-dx/2;
% tack on last bin end
x = [x x(end)+dx(end)];
% 
y = [y1 y1(:,end)]';
