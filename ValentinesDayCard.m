function ValentinesDayCard()

[x,y,z]=meshgrid(-3:0.03:3);
v=((x.^2)+(9/4).*(y.^2)+(z.^2)-1).^3-(x.^2).*(z.^3)-9/80.*(y.^2).*(z.^3);
isosurface(x,y,z,v,0);
rotate3d on
colormap(hsv(1));
daspect([.75 1 .75])
axis off
axis vis3d
camlight('headlight')
title('Happy Valentines Day Nati!')

end