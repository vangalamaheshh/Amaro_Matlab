function put_on_plane(h,p)

f=gcf;
a=gca;

ax = ancestor(h(1),'axes');
if isempty(ax) || ax==0,
  error(id('InvalidHandle'),'H must contain axes children only.');
end
ax=[get(ax,'xlim')' get(ax,'ylim')' [0 0]' ]; % get(ax,'zlim')'];
origin = sum(ax)/2;

%if diff(ax(:,3))~=0
%  error('should be a 2D plane');
%end

x=[ ax(1,:); ax(2,1) ax(1,2) ax(1,3); ax(2,1) ax(2,2) ax(1,3); ax(1,1) ax(2,2) ax(1,3);ax(1,:)];
xc=x-repmat(origin,5,1);

mid=[ mean(p(1:4,1)) mean(p(1:4,2)) mean(p(1:4,3))];
pc=p-repmat(mid,5,1);

pdx=pdist(x);
pdp=pdist(p);
scale=eye(3);
scale(1,1)=pdp(1)/pdx(1); % dist(1,2)
scale(2,2)=pdp(5)/pdx(5); % dist(2,3);

nx=null(xc(1:4,:));
if det([xc(1:2,:); nx']) < 0
  nx=-nx;
end

np=null(pc(1:4,:));
if det([pc(1:2,:); np']) < 0
  np=-np;
end

figure(f);
axes(a);
c=copyobj(h,a);
c(end+1)=patch(x(:,1),x(:,2),x(:,3),1);
%x
scaleshift(c,scale,[0 0 0],origin);
rp(:,1)=get(c(end),'Xdata');rp(:,2)=get(c(end),'Ydata');rp(:,3)=get(c(end),'Zdata');
%rp
xc=xc*scale;

% axis off;

if any(nx-np) % not same normal
  npx=null([nx np]');
  if det([nx'; np'; npx']) < 0
    npx=-npx;
  end
  
  nnpx=null(npx');
  x=[ nx np ]'*nnpx;
  y=x*x';
  angle=acos(y(1,2))*180/pi;
  rotate(c,npx',angle,origin);
end

rp(:,1)=get(c(end),'Xdata');rp(:,2)=get(c(end),'Ydata');rp(:,3)=get(c(end),'Zdata');
rp=rp-repmat(origin,size(rp,1),1);
nrp=null(rp(1:2,:));
if det([rp(1:2,:); nrp']) < 0
  nrp=-nrp;
end

x=[pc(2,:); rp(2,:)];
x=x./repmat(sqrt(sum(x.^2,2)),1,3);
y=x*x';
%nrp'
%np'
acos(y(1,2))*180/pi
rotate(c,np',acos(y(1,2))*180/pi,origin);

scale=1;

scaleshift(c,scale,mid-origin,origin);
delete(c(end));
%mid
% keyboard





