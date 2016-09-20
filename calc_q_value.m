function [q,pi0]=calc_q_value(pvec)

[sp,ord]=sort(pvec);
m=length(pvec);

x=0:0.01:0.95;
for i=1:length(x)
    y(i)=length(find(pvec>x(i)))/m/(1-x(i)+eps);
end 
%figure(1); clf;
%plot(x,y,'x'); hold on
%pp=spaps([x 1],[y y(end)],0.03,(1-[x 1]).^2,2);
%pnts=fnplt(pp,'r',[0 1],3); 
%pi0=pnts(end);


pol=polyfit(x,y,3);
pi0=min(polyval(pol,1),1);
if (0)
  figure(1); clf;
  plot(x,y,'x'); hold on
  plot(x,polyval(pol,x));
end
% keyboard
q(m)=pi0*sp(m);
for i=m-1:-1:1
    q(i)=min(pi0*m*sp(i)/i,q(i+1));
end 
q(ord)=q;
