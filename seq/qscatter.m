function qscatter(Ga,Gb)

X=[];
X.gene=Ga.gene;X.qa=Ga.q;X.qb=mapacross(X.gene,Gb.gene,Gb.q);

cap = 1e-6;
X.qa(X.qa<cap)=cap; X.qb(X.qb<cap)=cap; X.x=-log10(X.qa);X.y=-log10(X.qb);

x_sigthresh=1;
y_sigthresh=1;

clf,set(gca,'position',[0.1 0.14 0.8 0.8]),scatter(X.x,X.y,10,[0.8 0.8 0]),line(xlim,xlim,'color',[0 0 1])
line(xlim,[1 1]*y_sigthresh,'color',[1 0 0]);line([1 1]*x_sigthresh,ylim,'color',[1 0 0]);
gidx1=find(X.x>=x_sigthresh & X.y>=y_sigthresh); cols1 = repmat([0 0 0],length(gidx1),1);
gidx2=find(X.x<x_sigthresh & X.y>=y_sigthresh); cols2 = repmat([0 0.7 0],length(gidx2),1);
gidx3=find(X.x>=x_sigthresh & X.y<y_sigthresh); cols3 = repmat([1 0 0],length(gidx3),1);
gidx = [gidx1;gidx2;gidx3]; cols = [cols1;cols2;cols3];
textfit(X.x(gidx)+0.05,X.y(gidx),X.gene(gidx),'color',cols,'fontsize',6);
xlabel('old','fontsize',20);ylabel('new','fontsize',20);
