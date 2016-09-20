function ylabels_by_group(labels,varargin)
% ylabels_by_group(labels)
%
% Given an array of text labels,
% finds the boundaries between identical labels,
% and draws the ticks there, and the labels in between.
%
% Mike Lawrence 2009-07-17

if isnumeric(labels), labels = num2cellstr(labels); end

l1 = [labels(1:end)];
l2 = [labels(2:end);'***'];
b = find(~strcmp(l1,l2));

set(gca,'ytick',[0;b+0.5],'yticklabel',{});

xl = xlim;
xpos = xl(1) - 0.03*(xl(2)-xl(1));

for i=1:length(b)
  if i==1, y1=0; else y1=b(i-1); end
  y2 = b(i);
  ypos = y1 + 0.1*(y2-y1);
  text(xpos,ypos,labels{b(i)},'clipping','off',...
    'verticalalign','top',varargin{:},'horizontalalign','right','interpreter','none');
end

