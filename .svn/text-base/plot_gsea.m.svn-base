function plot_gsea(ro_set,n)

plot((0:n)/n,(0:n)/n,'b-');
hold on;
vec=zeros(n,1);
vec(ro_set)=1;
cs=cumsum([0; vec]);
plot((0:n)/n,cs/cs(end),'r-');
cs(end)
