function draw_context_plot(R,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'fontsize',20);

R = reorder_struct(R,grep('total',R.regulatory_potential,1));
R = reorder_struct(R,grep('total',R.transcript_zone,1));

close,figure(1),clf
barweb(R.mutrate',R.ci_low',R.ci_high',0.8,{'mutation type'});
ylabel('mutations per million sites','fontsize',P.fontsize); title('Mutation rates','fontsize',P.fontsize);
legend(R.base_context,'location','eastoutside','interpreter','none');
set(gca,'fontsize',P.fontsize);
set(gcf,'color',[1 1 1]);
ymax = max([R.mutrate(:);R.ci_low(:);R.ci_high(:)]);
ylim([0 ymax*1.05]);
