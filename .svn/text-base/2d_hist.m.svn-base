function hd=2d_hist(x,y,rx,ry)

xn=zeros(size(x));
yn=zeros(size(x));
for i=1:(length(rx)-1)
  xn(find( x >=rx(i) & x <rx(i+1)))=i;
  yn(find( y >=ry(i) & y <ry(i+1)))=i;
end

hd=crosstab_full(xn,yn,1:max(length(rx),length(ry)));
