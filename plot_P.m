function [lh,ph]=plot_P(P,n,f)

sP=sort(P);
x=sP(1:n);
for i=1:n;
  s=sum(P(:,2:end)<=sP(i,1));
  [r1(i),r2(i)]=two_sided_prctile(s,f);
end

clf;
keyboard
ph=patch([ x fliplr(x) x(1)],[r1 fliplr(r2) r1(1)],[0.7 0.7 0.7]);
set(ph,'LineStyle','none','FaceAlpha',0.5);
hold on
lh=plot(x,1:n);
