function put_in_3d(h,p,flip_ax,fig_off)

if ~exist('fig_off','var')
  fig_off=0;
end

f=gcf;
a=gca;

axs = ancestor(h(1),'axes');
if isempty(axs) || axs==0,
  error(id('InvalidHandle'),'H must contain axes children only.');
end
ax=[get(axs,'xlim')' get(axs,'ylim')' get(axs,'zlim')'];
origin = sum(ax)/2;

v=get(axs,'View');
if (v(1)~=0 && v(2)~=90)
  origin(3)=ax(1,3);
end;

szx=diff(ax);

mirror=0;
if size(p,1)==5 % plane
  ax(:,3)=origin(3);
  pd=dist(p,[],'euclid');
  szp=[ pd(1,2) pd(2,3) 0];
  p4=p(1:4,:);
  mid=mean(p4,1);
elseif size(p,1)==3 % plane basis
  ax(:,3)=origin(3);
  pd=dist(p,[],'euclid');
  szp=[ pd(1,2) pd(1,3) 0]; 
  mid=mean(p(2:3,:),1);
  p=[ p(1,:); p(2,:); p(2,:)+p(3,:)-p(1,:); p(3,:); p(1,:)];
elseif size(p,1)==4 % 3d basis
  pd=dist(p,[],'euclid');
  szp=[ pd(1,2) pd(1,3) pd(1,4)];   
  mid=mean(p(2:3,:),1);
  p=[ p(1,:); p(2,:); p(2,:)+p(3,:)-p(1,:); p(3,:); p(1,:)];
else
  error('p has too few or too many rows');
end


  
%if diff(ax(:,3))~=0
%  error('should be a 2D plane');
%end

x=[ ax(1,:); ax(2,1) ax(1,2) ax(1,3); ax(2,1) ax(2,2) ax(1,3); ax(1,1) ax(2,2) ax(1,3);ax(1,:)];
xc=x-repmat(origin,5,1);

pc=p-repmat(mid,5,1);

scale=zeros(3,3);
scale=diag(szp./szx);
if szp(3)==0
  scale(3,3)=1;
end
% diag(scale)

if exist('flip_ax','var')
  scale=scale.*(ones(3,3)-2*diag(flip_ax));
%   scale
end


nx=null(xc(1:4,:));
if det([xc(1:2,:); nx']) < 0
  nx=-nx;
end

np=null(pc(1:4,:));
if det([pc(1:2,:); np']) < 0
  np=-np;
end

figure(f);
if fig_off
  set(gcf,'Visible','off');
end
axes(a);
if fig_off
  set(gcf,'Visible','off');
end
c=copyobj(h,a);
if fig_off
  set(gcf,'Visible','off');
end
c(end+1)=patch(x(:,1),x(:,2),x(:,3),1);
if fig_off
  set(gcf,'Visible','off');
end
%x
scaleshift(c,scale,[0 0 0],origin);
rp(:,1)=get(c(end),'Xdata');rp(:,2)=get(c(end),'Ydata');rp(:,3)=get(c(end),'Zdata');
%rp
xc=xc*scale;

nx=null(xc(1:4,:));
if det([xc(1:2,:); nx']) < 0
  nx=-nx;
end

% axis off;

if any(nx-np) % not same normal
  if ~any(nx+np) % nx=-np;
    rotate(c,sum(xc(1:2,:),1),180,origin);
  else
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
y(y>1)=1;
% acos(y(1,2))*180/pi
% y(1,2)
rotate(c,np',acos(y(1,2))*180/pi,origin);

scale=1;

scaleshift(c,scale,mid-origin,origin);
delete(c(end));
%mid
% keyboard





