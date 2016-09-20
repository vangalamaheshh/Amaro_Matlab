function ct=generate_mesh(x,y,xvals,yvals)

xi=x;
for i=1:(length(xvals)-1)
  xi(x>=xvals(i) & x<xvals(i+1))=i;
end

yi=y;
for i=1:(length(yvals)-1)
  yi(y>=yvals(i) & y<yvals(i+1))=i;
end

ct=crosstab_full(xi,yi,1:(length(xvals)-1));
