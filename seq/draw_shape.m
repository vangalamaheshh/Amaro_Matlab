function draw_shape(shp,x,y,col,style,halo_col)
%
% private subfunction of draw_genefig
%

% shape definitions

if strcmp(shp,'sqr')
  shp = struct('x',[0 1 1 -1 -1 0],'y',[-1 -1 1 1 -1 -1],'s1',0.013,'s2',0.025,'s3',0.04);
elseif strcmp(shp,'cir')
  shp = struct('x',cos(0:pi/12:2*pi),'y',sin(0:pi/12:2*pi),'s1',0.015,'s2',0.028,'s3',0.045);
elseif strcmp(shp,'dmd')
  shp = struct('x',[0 0.9 0 -0.9 0],'y',[1.2 0 -1.2 0 1.2],'s1',0.016,'s2',0.036,'s3',0.05);
elseif strcmp(shp,'tri')
  shp = struct('x',[0 1.2 0 -1.2 0],'y',[-0.7 -0.7 1.3 -0.7 -0.7],'s1',0.015,'s2',0.031,'s3',0.047);
elseif strcmp(shp,'inv')
  shp = struct('x',[0 1.2 0 -1.2 0],'y',[0.7 0.7 -1.3 0.7 0.7],'s1',0.015,'s2',0.031,'s3',0.047);
elseif strcmp(shp,'dot')
  shp = struct('x',cos(0:pi/4:2*pi),'y',sin(0:pi/4:2*pi),'s1',0.005,'s2',0.008,'s3',0.02);
else
  error('unknown shape');
end

if exist('style','var')
  if strcmp(style,'small')
    patch(x+shp.x*(shp.s2/2),y+shp.y*(shp.s2/2), halo_col, 'linestyle','none');
  elseif strcmp(style,'large')
     patch(x+shp.x*(shp.s3/2),y+shp.y*(shp.s3/2), halo_col, 'linestyle','none');
  elseif strcmp(style,'outlined')
    patch(x+shp.x*((shp.s1+0.008)/2),y+shp.y*((shp.s1+0.008)/2), halo_col, 'linestyle','none');
  end
end

patch(x+shp.x*(shp.s1/2),y+shp.y*(shp.s1/2), col, 'linestyle','none');

end
