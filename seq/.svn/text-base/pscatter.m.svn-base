function pscatter(Ga,Gb)

if nargin==1
  if isfield(Ga,'pa') && isfield(Ga,'pb')
    Gb=Ga;
    Ga.p = Ga.pa;
    Gb.p = Gb.pb;
  elseif isfield(Ga,'p1') && isfield(Ga,'p2')
    Gb=Ga;
    Ga.p = Ga.p1;
    Gb.p = Gb.p2;
  else
    error('needs two p-values');
  end
end

X=[];X.gene=Ga.gene;X.pa=Ga.p;X.pb=mapacross(X.gene,Gb.gene,Gb.p);
X.pa(X.pa==0)=1e-15; X.pb(X.pb==0)=1e-15; X.x=-log10(X.pa);X.y=-log10(X.pb);
X.qa=calc_fdr_value(X.pa); [tmp ord] = sort(abs(X.qa-0.1)); x_sigthresh=X.x(ord(1));
X.qb=calc_fdr_value(X.pb); [tmp ord] = sort(abs(X.qb-0.1)); y_sigthresh=X.y(ord(1));
clf,set(gca,'position',[0.1 0.14 0.8 0.8]),scatter(X.x,X.y,10,[0.8 0.8 0]),line(xlim,xlim,'color',[0 0 1])
line(xlim,[1 1]*y_sigthresh,'color',[1 0 0]);line([1 1]*x_sigthresh,ylim,'color',[1 0 0]);
gidx1=find(X.x>=x_sigthresh & X.y>=y_sigthresh); cols1 = repmat([0 0 0],length(gidx1),1);
gidx2=find(X.x<x_sigthresh & X.y>=y_sigthresh); cols2 = repmat([0 0.7 0],length(gidx2),1);
gidx3=find(X.x>=x_sigthresh & X.y<y_sigthresh); cols3 = repmat([1 0 0],length(gidx3),1);
gidx = [gidx1;gidx2;gidx3]; cols = [cols1;cols2;cols3];
textfit(X.x(gidx)+0.05,X.y(gidx),X.gene(gidx),'color',cols,'fontsize',8);
xlabel('old','fontsize',20);ylabel('new','fontsize',20);
