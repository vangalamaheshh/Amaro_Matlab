function plot_perlane_chunk_correlations(X)
% Mike Lawrence 2009-07-17

X.corr = 1 - dist(X.relcap',[],'correlation');
clf;imagesc(X.corr);colormap(jet);colorbar
xlabels_by_group(X.lane.short,'fontsize',7);
ylabels_by_group(X.lane.short,'fontsize',7);
tnline = 0.5+find(strcmp('N',X.lane.tn(1:end-1)) & strcmp('T',X.lane.tn(2:end)));
if length(tnline)==1
  line([tnline tnline],[0 X.nlanes+0.5],'color',[0 0 0],'linewidth',2);
  line([0 X.nlanes+0.5],[tnline tnline],'color',[0 0 0],'linewidth',2);
end
set(gca,'position',[0.05 0.05 0.80 0.85])
