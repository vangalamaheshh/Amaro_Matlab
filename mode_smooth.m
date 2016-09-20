function [m,sg]=mode_smooth(x,sig,nb)

hx=(min(x)-3*sig):(range(x)+6*sig)/(nb-1):(max(x)+3*sig);

g=zeros(length(x),length(hx));
for i=1:length(x)
  g(i,:)=normpdf(hx,x(i),sig);
end
sg=sum(g);
[msg,msgi]=max(sg);
m=hx(msgi);
