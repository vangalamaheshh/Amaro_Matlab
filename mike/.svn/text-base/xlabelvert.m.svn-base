function hText = xlabelvert(varargin)

XTick = get(gca,'xtick');
if length(XTick)==1, fprintf('note: bug in xticklabel for length(XTick)==1\n'); return; end

XTickLabel = get(gca,'xticklabel');
rot = 90;


xticklabel_rotate(XTick,rot,XTickLabel,'interpreter','none',varargin{:})   

if nargout==0
  clear hText
end
