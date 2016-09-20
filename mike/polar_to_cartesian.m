function [x y] = polar_to_cartesian(r,theta)

x = r.*cos(theta);
y = r.*sin(theta);

