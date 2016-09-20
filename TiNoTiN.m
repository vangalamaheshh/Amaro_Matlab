tau=0:0.1:15;  % observed total copy ratio range from ACS 
TiN=0.5;         % some TiN 
ft=0.05;            % some hefty tumor het allele shift from 0.5
Ta=ft*tau;
TiNo=1/ft*((TiN*Ta+1-TiN)./(2-2*TiN-tau*TiN));
plot(tau,TiNo)
xlabel('tau')
ylabel('TiNo')


K=(2/TiN)-2;
y=tau./(K+tau);
plot(tau,y)
xlabel('tau')
ylabel('MM TiNo')

b.theoretical_TiN=2./((b.tau./b.modal_TiN)-b.tau+2);

cfilt=reorder_struct(call_filtered,call_filtered.x>=b.xs(seg)&call_filtered.x<=b.xe(seg));
figure()
plot(b.tau,2./((b.tau./b.modal_TiN)-b.tau+2),'b.','MarkerSize',20)

figure()
subplot(2,1,1)
plot(cfilt.normal_f,cfilt.tumor_f,'b.','MarkerSize',20)
xlabel('Normal F','FontSize',20)
ylabel('Tumor F','FontSize',20)
subplot(2,1,2)
plot(cfilt.x,cfilt.normal_f,'r.','MarkerSize',20)
hold on
plot(cfilt.x,cfilt.tumor_f,'b.','MarkerSize',20)
ylabel('AF','FontSize',20)
xlabel('Position','FontSize',20)