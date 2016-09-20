function draw_CpG_island_rate_plot(R,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'fontsize',20);

R = reorder_struct(R,grep('IGR|intron|exon',R.transcript_zone,1));
Rtot = reorder_struct(R,grep('total',R.base_context,1));
Rc = reorder_struct(R,grep('CpG_transition',R.base_context,1));
Rnc = Rtot;
Rnc.base_context = repmat({'non_CpG_transition'},slength(Rnc),1);
Rnc.n_muts = Rtot.n_muts - Rc.n_muts;
[rate confint] = binofit(Rnc.n_muts,Rnc.Ncov);
Rnc.mutrate = 1e6*rate;
Rnc.ci_low = 1e6*confint(:,1);
Rnc.ci_high = 1e6*confint(:,2);
Rnc.stdev = (Rnc.ci_high-Rnc.mutrate)/1.96;

xlabels = {'island','shore','sea','total'};

close,figure(1),clf
subplot('position',[0.3 0.55 0.63 0.4]);
barweb(reshape(Rtot.mutrate,3,4)',reshape(Rtot.ci_low,3,4)',reshape(Rtot.ci_high,3,4)',0.8,xlabels);
ylabel('mutations per million sites','fontsize',P.fontsize); title('All mutations','fontsize',P.fontsize);
legend({'IGR','intron','exon'},'location','eastoutside');
set(gca,'fontsize',P.fontsize);
subplot('position',[0.1 0.05 0.4 0.4]);
barweb(reshape(Rc.mutrate,3,4)',reshape(Rc.ci_low,3,4)',reshape(Rc.ci_high,3,4)',0.8,xlabels);
ylabel('mutations per million sites','fontsize',P.fontsize); title('CpG transitions','fontsize',P.fontsize);
set(gca,'fontsize',P.fontsize);
subplot('position',[0.58 0.05 0.4 0.4]);
barweb(reshape(Rnc.mutrate,3,4)',reshape(Rnc.ci_low,3,4)',reshape(Rnc.ci_high,3,4)',0.8,xlabels);
ylabel('mutations per million sites','fontsize',P.fontsize); title('mutations other than CpG transitions','fontsize',P.fontsize);
set(gca,'fontsize',P.fontsize);
set(gcf,'position',[78 94 1154 851],'color',[1 1 1]);

