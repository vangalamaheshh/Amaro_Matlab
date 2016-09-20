N=500;
yt1=betarnd(96,5,N,1);
yt2=betarnd(5,96,N,1);
yn1=betarnd(61,41,N,1);
yn2=betarnd(41,61,N,1);
x=1:N;
b=0:0.01:1
nt1=hist(yt1,b);
nt2=hist(yt2,b);
nn1=hist(yn1,b);
nn2=hist(yn2,b);

subplot(1,10,1:8)
plot(x,yt1,'b.',x,yt2,'b.',x,yn1,'r.',x,yn2,'r.')
subplot(1,10,9:10)
stairs(b,nt1); hold on; stairs(b,nt2); stairs(b,nn1,'color','r'); stairs(b,nn2,'color','r'); hold off
view(90,90)
set(gca,'xticklab',[])