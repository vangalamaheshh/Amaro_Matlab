function heatmap(y,x,numybins,numxbins,logflag)
% heatmap(y,x,numybins,numxbins,logflag)
% Mike Lawrence 2009-2010

if ~exist('numxbins','var'), numxbins = 20; end
if ~exist('numybins','var'), numybins = 20; end
if ~exist('logflag','var'), logflag = false; end

if length(x) ~= length(y), error('x and y must be same length'); end

Z = zeros(numybins+1,numxbins+1);
xmin = min(x); xmax = max(x); xspan=xmax-xmin; xbinsize = xspan/numxbins;
ymin = min(y); ymax = max(y); yspan=ymax-ymin; ybinsize = yspan/numybins;

binx = 1+round((x-xmin)/xbinsize);
biny = 1+round((y-ymin)/ybinsize);

for i=1:length(x), Z(biny(i),binx(i)) = Z(biny(i),binx(i)) + 1; end

if logflag, Z = log10(Z); end

imagesc(Z);
colorbar;
l={}; for i=xmin+xbinsize/2:xbinsize:xmax, l=[l;num2str(round(i))]; end
set(gca,'xtick',0.5+[1:numxbins],'xticklabel',l);

l={}; for i=ymin+ybinsize/2:ybinsize:ymax, l=[l;num2str(round(i))]; end
set(gca,'ytick',0.5+[1:numybins],'yticklabel',l);

if ~logflag
  for x=1:numxbins+1, for y=1:numybins+1
      if Z(y,x)>0 & Z(y,x)<10, text(x,y,num2str(Z(y,x))); end
  end,end
end

