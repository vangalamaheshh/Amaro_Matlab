nt=10000;
N=10;
pv=zeros(nt,1);
pv2=zeros(nt,1);
md=zeros(nt,1);
for i=1:nt
    r=randperm(N);
    pv(i)=ranksum_many(1:N,r(1:(N/2)),r(((N/2)+1):end),-1);
    pv2(i)=ttest2_many(1:N,r(1:(N/2)),r(((N/2)+1):end),-1);
    md(i)=mean(r(1:(N/2)))-mean(r(((N/2)+1):end));
end


figure(1); clf;
subplot(2,1,1);
hist(pv,0:0.01:1);
subplot(2,1,2);
hist(pv2,0:0.01:1);

figure(2); clf;
plot(md,pv,'x');
