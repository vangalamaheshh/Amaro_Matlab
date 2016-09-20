function scaleshift(h,scale,shiftvec,origin)

if nargin < 4
  ax = ancestor(h(1),'axes');
  if isempty(ax) || ax==0,
    error(id('InvalidHandle'),'H must contain axes children only.');
  end
  origin = sum([get(ax,'xlim')' get(ax,'ylim')' get(ax,'zlim')'])/2;
end

if numel(shiftvec)==2
  shiftvec(3)=0;
end

for i=1:numel(h),
  t = get(h(i),'type');
  skip = 0;
  if strcmp(t,'surface') || strcmp(t,'line') || strcmp(t,'patch')
    
    % If patch, rotate vertices  
    if strcmp(t,'patch')
       verts = get(h(i),'Vertices');
       x = verts(:,1); y = verts(:,2); 
       if size(verts,2)>2
          z = verts(:,3);
       else
          z = [];
       end
       
    % If surface or line, rotate {x,y,z}data   
    else
       x = get(h(i),'xdata');
       y = get(h(i),'ydata');
       z = get(h(i),'zdata');
    end
    
    if isempty(z)
       z = -origin(3)*ones(size(y));
    end
    [m,n] = size(z);
    if numel(x) < m*n
      [x,y] = meshgrid(x,y);
    end
  elseif strcmp(t,'text')
    p = get(h(i),'position');
    x = p(1); y = p(2); z = p(3);
  elseif strcmp(t,'image')
    x = get(h(i),'xdata');
    y = get(h(i),'ydata');
    z = zeros(size(x));
  else
    skip = 1;
  end
  
  if ~skip,
    [m,n] = size(x);
    
    newxyz = [x(:)-origin(1),y(:)-origin(2),z(:)-origin(3)];
    newxyz = newxyz*scale;
    newx = origin(1)+shiftvec(1)+reshape(newxyz(:,1),m,n);
    newy = origin(2)+shiftvec(2)+reshape(newxyz(:,2),m,n);
    newz = origin(3)+shiftvec(3)+reshape(newxyz(:,3),m,n);
    if strcmp(t,'surface') || strcmp(t,'line')
      set(h(i),'xdata',newx,'ydata',newy,'zdata',newz);
    elseif strcmp(t,'patch')
      set(h(i),'Vertices',[newx,newy,newz]); 
    elseif strcmp(t,'text')
      set(h(i),'position',[newx newy newz])
    elseif strcmp(t,'image')
      set(h(i),'xdata',newx,'ydata',newy)
    end
  end
end

function str=id(str)
str = [':shift:' str];

