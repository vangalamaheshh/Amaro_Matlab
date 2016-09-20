function y=min_over_range(x,dim,min_val,max_val)

x(x<min_val)=NaN;
x(x>max_val)=NaN;
y=nanmin(x,[],dim);
