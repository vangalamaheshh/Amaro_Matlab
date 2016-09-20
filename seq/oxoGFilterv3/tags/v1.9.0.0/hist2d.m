function [counto,xco,yco,ho] = hist2d(x,y,nx,ny,CLIM1,P)
% HIST2D calculates a 2-dimensional histogram
%    [N,xb,yb,h] = HIST2D(X,Y,NX,NY,CLIM,P) bins the X and Y data into bins NX and NY in
%    both dimensions. The color (z) scale is set by CLIM (optional) and a
%    colorbar is displayed unless P.BAR==false (optional).  counto is the count in each bin.  
%
%    N = HIST2D(X,Y,M), where M is a scalar, bins the data into M equally
%    spaced bins in both dimensions
%
%    N = HIST2D(X,Y,B), where B is a vector, bins the data with centers
%    specified by B 
%
%    The number of bins or centers can be specified individually for either
%    dimension using N = HIST2D(X,Y,NX,NY) or N = HIST2D(X,Y,BX,BY)
%
%    [N,BX,BY] = HIST2D(...) also returns the positions of the bin centers
%    as two matrices in BX and BY
%
%    HIST2D(...) without output arguments produces a colormapped image plot
%    of the 2d histogram
%
%    HIST2D(...) also returns the positions of the bin centers
%    as two matrices in BX and BY
%
%    The last input argument is a struct of options. 
%       P.BAR=true (default) will make the colorbar 
%       P.MARGINHISTS=true (default)will make the margin histos.
%
%    The last output argument is a vector of handles to
%      ho(1) = 2d imagesc
%      ho(2) = colorbar
%      ho(3) = vertical histo on left
%      ho(4) = horizontal histo on bottom
%
% EXAMPLE
%   yin = randn(1,1000);
%   xin = randn(1,1000);
%   [n,x,y] = hist2d(xin,yin,11);
%   imagesc(x(1,:),y(:,1),n); hold on; plot(xin,yin,'y.'); colorbar

if ~exist('nx')
   nx = 10;
end

if ~exist('ny')
   ny = nx;
end

if ~exist('CLIM1')
   CLIM1 = NaN*[1 100];
end

if nargin<6
    P.BAR=true;
    P.MARGINHISTS=true;
end

if ~isfield(P,'BAR')
    P.BAR=true;
end    
if ~isfield(P,'MARGINHISTS')
    P.MARGINHISTS=true;
end

if length(x) ~= length(y)
   error(sprintf('x and y must be same size ( %g ~= %g )',length(x),length(y)));
end
   
[dummy,xc] = hist(x,nx);
[dummy,yc] = hist(y,ny);

count = [];

for i = 1:length(yc)
   if i == 1
      lbound = -Inf;
   else
      lbound = (yc(i-1) + yc(i)) / 2;
   end
   if i == length(yc)
      ubound = inf;
   else
      ubound = (yc(i) + yc(i+1)) /2;
   end
   count(i,:) = hist(x((y >= lbound) & (y < ubound)),xc);
end

[xc, yc] = meshgrid(xc, yc);

hh=NaN;

if (nargout == 0)|| (nargout>3)
    %cmap=1-(1-jet).*((gray).^0.25);
    v=(0:254)'/254;
    cmap=[v 0*v 1-v]; cmap=[[1 1 1]; cmap];
    zc=count; zc(zc<1)=0.1;
    if isnan(CLIM1(1))
        CLIM1=[0.9 10^(round(log10(max(zc(:))))+1)];
    end
    if (CLIM1(1)==1)
        CLIM1(1)=0.9;
    end
    if (CLIM1(2)==1)
        CLIM1(2)=2;
    end
    if (P.MARGINHISTS)
        subplot('position',[0.25 0.28 0.65 0.65])
    end
    h1=imagesc(xc(1,:),yc(:,1),log10(zc),log10(CLIM1));
    hh(1)=get(h1,'parent');
    axis xy; colormap(cmap); grid on;
    if (P.BAR)
        hh(2)=colorbar;
        zt=get(hh(2),'ytick');
        lzt=round(10.^zt);
        zt=log10(unique(lzt));
        set(hh(2),'ytick',zt,'yticklabel',num2str(round(10.^zt')));
    end
    if (P.MARGINHISTS)
        set(hh(1),'xticklabel',[],'yticklabel',[])
        dx=diff(xc(1,1:2)); xlim1=[xc(1,1)-dx/2 xc(1,end)+dx/2];
        dy=diff(yc(1:2,1)); ylim1=[yc(1,1)-dy/2 yc(end,1)+dy/2];
        subplot('position',[0.11 0.28 0.12 0.65])
        h1(3)=barh(yc(:,1),sum(count,2),'hist');
        hh(3)=get(h1(3),'parent');
        ylim(ylim1)
        p1=get(hh(1),'position');
        subplot('position',[0.25 0.117 p1(3) 0.12])
        h1(4)=bar(xc(1,:),sum(count,1),'hist');
        hh(4)=get(h1(4),'parent');
        xlim(xlim1)
        set(h1(3:4),'FaceColor',0.5*[1 1 1])
    end
end

if (nargout>0)
    counto = count;
    xco = xc;
    yco = yc;
    ho=hh;
end

%%
function test
x=randn(5000,1);
y=randn(5000,1);
hist2d(x,y,-3:0.2:3,-3:0.2:3)
