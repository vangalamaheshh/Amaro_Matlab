p1=[0 0 0; 2 0 0; 2 1 0; 0 1 0; 0 0 0];
R=[10 50 2; 0 7 10; 0 0 19];
R=R+R';
[V,D]=eig(R);
q=[ 3 5 2];
mid1=mean(p1(1:4,:));
p2=(p1-repmat(mid1,size(p1,1),1))*V+repmat(q+mid1,size(p1,1),1);

figure(1); clf;
ph(1)=patch(p1(:,1),p1(:,2),p1(:,3),1);
ph(2)=patch(p2(:,1),p2(:,2),p2(:,3),2);
set(ph,'FaceAlpha',0.2);
grid on

mid1=mean(p1(1:4,:));
mid2=mean(p2(1:4,:));

p2=p2-repmat(mid2,size(p2,1),1);
p1=p1-repmat(mid1,size(p1,1),1);

figure(1); clf;
ph(1)=patch(p1(:,1),p1(:,2),p1(:,3),1);
ph(2)=patch(p2(:,1),p2(:,2),p2(:,3),2);
set(ph,'FaceAlpha',0.2);
grid on

n1=null(p1(1:4,:));
if det([p1(1:2,:); n1']) < 0
  n2=-n2;
end
n2=null(p2(1:4,:));
if det([p2(1:2,:); n2']) < 0
  n2=-n2;
end

line([ 0 n1(1)],[0 n1(2)],[0 n1(3)]); hold on;
line([ 0 n2(1)],[0 n2(2)],[0 n2(3)]);
n12=null([n1 n2]');
if det([n1'; n2'; n12']) < 0
  n12=-n12;
end

nn12=null(n12');
line([0 n12(1)],[0 n12(2)],[0 n12(3)],'Color','r');
x=[ n1 n2 ]'*nn12;
y=x*x';
a=acos(y(1,2))*180/pi;
rotate(ph(2),n12',-a,[ 0 0 0]);
rp(:,1)=get(ph(2),'Xdata');rp(:,2)=get(ph(2),'Ydata');rp(:,3)=get(ph(2),'Zdata');
x=[p1(2,:); rp(2,:)];
line([ 0 x(1,1)],[0 x(1,2)],[0 x(1,3)]); hold on;
line([ 0 x(2,1)],[0 x(2,2)],[0 x(2,3)]); hold on;

x=x./repmat(sqrt(sum(x.^2,2)),1,3);
y=x*x';
rotate(ph(2),n1',-acos(y(1,2))*180/pi,[ 0 0 0]);
