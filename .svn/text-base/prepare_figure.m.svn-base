function prepare_figure(fname,columns)

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 5 7.6 7.6*3/4]);
set(gcf,'Position',[520 690 560 420]);

hgsave(gcf,[ fname '.fig']);
print('-depsc',['-f' num2str(gcf)],[ fname '.eps']);
print('-dpng',['-f' num2str(gcf)],'-r90',[ fname '.png']);

