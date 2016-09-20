function sz=calc_point_size

pp=get(gcf,'PaperPosition');
wh=pp(3:4);

ap=get(gca,'Position');
xw=range(get(gca,'XLim'));
yw=range(get(gca,'YLim'));

awh=ap(3:4);
sz=[xw yw]./(awh.*wh*72);

