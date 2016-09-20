function [a b r] = linreg(x,y)
% Mike Lawrence 2009-10-02

x = as_row(x);
y = as_row(y);

good = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
n = sum(good,2);

x(~good) = 0;
y(~good) = 0;

xx = x.^2;
yy = y.^2;
xy = x.*y;

sx = sum(x,2);
sy = sum(y,2);
sxx = sum(xx,2);
syy = sum(yy,2);
sxy = sum(xy,2);

ntor = (n.*sxy)-(sx.*sy);
r = ntor ./ (sqrt((n.*sxx)-(sx.^2)) .* sqrt((n.*syy)-(sy.^2)));
a = ntor ./ ((n.*sxx)-(sx.^2));
b = (sy-(a.*sx)) ./ n;
