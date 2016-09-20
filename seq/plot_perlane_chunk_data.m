function plot_perlane_chunk_data(X)
% Mike Lawrence 2009-07-17

figure(1);clf;
if isfield(X.lane,'cov')
  subplot('position',[0.05 0.2 0.9 0.73]);
end
imagesc(X.relcap);
ylabels_by_group(chrlabel(X.chunk.chr));
xlabels_by_group(X.lane.short,'fontsize',7);
tnline = 0.5+find(strcmp('N',X.lane.tn(1:end-1)) & strcmp('T',X.lane.tn(2:end)));
if length(tnline)==1
  line([tnline tnline],[0 X.nchunks+0.5],'color',[0 0 0],'linewidth',2);
  text(tnline/2,X.nchunks*1.02,'normals','horizontalalign','center','fontsize',15);
  text((tnline+X.nlanes)/2,X.nchunks*1.02,'tumors','horizontalalign','center','fontsize',15);
end
colormap(bwr);colorbar;p1 = get(gca,'position');

if isfield(X.lane,'cov')
  subplot('position',[0.05 0.05 0.9 0.12]);bar(X.lane.cov/1e6,0.5);ylabel('reads x 10^6');
  xlabels_by_group(X.lane.short,'fontsize',7,'rotation',-90,'verticalalign','bottom');
  p2 = get(gca,'position');p2([1 3]) = p1([1 3]);set(gca,'position',p2);
  xlim([0.5 X.nlanes+0.5]);
end
