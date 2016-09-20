function xlabels_by_group(labels,lines_color,varargin)
% xlabels_by_group(labels,lines_flag,varargin)
%
% Given an array of text labels,
% finds the boundaries between identical labels,
% and draws the ticks there, and the labels in between.
%
% Mike Lawrence 2009-07-17

if ~exist('lines_color','var'), lines_color = []; end

if isnumeric(labels), labels = num2cellstr(labels); end

l1 = [labels(1:end)];
l2 = [labels(2:end);'***'];
b = find(~strcmp(l1,l2));

set(gca,'xtick',[0;b+0.5],'xticklabel',{});

yl = ylim();
if strcmp(get(gca,'ydir'),'reverse')
  ypos = yl(2);
else
  ypos = yl(1);
end

for i=1:length(b)
  if i==1
    xpos = (0 + b(i)) / 2;
  else
    xpos = (b(i-1) + b(i)) / 2;
  end
  text(xpos,ypos,labels{b(i)},'rotation',90,'clipping','off',...
    'horizontalalignment','right','interpreter','none',varargin{:});
end

if ~isempty(lines_color)
  for i=1:length(b)-1
    line([b(i) b(i)]+0.5,ylim,'color',lines_color);
  end
end

