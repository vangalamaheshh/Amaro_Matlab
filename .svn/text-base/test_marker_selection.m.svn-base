Ng=5000;
Ns0=5;
Ns1=20;
x=-1;

dat=[random('lognormal',0,1,Ng,Ns0) ...
     random('lognormal',x,1,Ng,Ns1)];

hist(dat(:),50) 
hist(log2(dat(:)),50) 


close all
%figure
%imagesc(dna_norm(dat)) 

[P,rs,P2]=marker_selection(dat,1:Ns0,Ns0+(1:Ns1),'snr',1000);


hist(P,20)

cls0=1:Ns0;
cls1=Ns0+(1:Ns1);
m0=mean(dat(:,cls0),2);
m1=mean(dat(:,cls1),2);
s0=std(dat(:,cls0),0,2);
s1=std(dat(:,cls1),0,2);


figure(1); clf;
subplot(2,1,1);
hist(P,20);
subplot(2,1,2);
hist(P2,20);


xs=[-1.5:0.5:1.5];
bs=zeros(length(xs),1);
bs2=zeros(length(xs),1);
Ps={};
inter=[];
for i=1:length(xs)
  Ng=10000;
%  Ns0=5;
%  Ns1=50;
   Ns0=20;
   Ns1=20;
  x=xs(i);
  
  dat=[random('lognormal',0,1,Ng,Ns0) ...
       random('lognormal',x,1,Ng,Ns1)];
%  dat=log2(dat);
  [P,rs,P2]=marker_selection(dat,1:Ns0,Ns0+(1:Ns1),'snr',1000);
  Ps{i,1}=P;
  Ps{i,2}=P2;
  bs(i)=length(find(P<0.05));
  bs2(i)=length(find(P2<0.05));
  inter(i)=length(intersect(find(r2.Ps{i,1}<0.05),...
                            find(r2.Ps{i,2}<0.05)));

  disp(xs(i));
  figure(1); clf;
  plot(xs(1:i),bs(1:i)); hold on
  plot(xs(1:i),bs2(1:i),'r-');  
end
save power_test_lognormal_balanced.mat Ps bs bs2 xs inter 


r1=load('power_test_lognormal.mat');
r2=load('power_test_lognormal1.mat');

figure(1); clf;
plot(r1.xs,r1.bs,'b-o'); hold on
plot(r1.xs,r1.bs2,'r-o');  
legend({'2*min(p,1-p)','1-sided abs(score)'});
plot(r2.xs,r2.bs,'b-o'); hold on
plot(r2.xs,r2.bs2,'r-o');  

fs=18;
fs1=14;
title('log normal, snr statistic, n0=5, n1=50','FontSize',fs);
xlabel('\Delta','FontSize',fs);
ylabel('Power','FontSize',fs);
set(gca,'FontSize',fs1);

print -dpdf -f1 lognormal_power.pdf

inter=[];
for i=1:length(r1.xs)
  inter(i)=length(intersect(find(r1.Ps{i,1}<0.05),find(r1.Ps{i,2}< ...
                                                    0.05)));
end
r1.inter=inter;

inter=[];
for i=1:length(r2.xs)
  inter(i)=length(intersect(find(r2.Ps{i,1}<0.05),find(r2.Ps{i,2}< ...
                                                    0.05)));
end
r2.inter=inter;

[ r1.bs'; r1.bs2'; r1.inter]

save power_test_lognormal_runs.mat r1 r2

r1=load('power_test_normal.mat');
figure(1); clf;
plot(r1.xs,r1.bs,'b-o'); hold on
plot(r1.xs,r1.bs2,'r-o');  
legend({'2*min(p,1-p)','1-sided abs(score)'});
%plot(r2.xs,r2.bs,'b-o'); hold on
%plot(r2.xs,r2.bs2,'r-o');  

fs=18;
fs1=14;
title('normal, snr statistic, n0=5, n1=50','FontSize',fs);
xlabel('\Delta','FontSize',fs);
ylabel('Power','FontSize',fs);
set(gca,'FontSize',fs1);

print -dpdf -f1 normal_power.pdf


r1=load('power_test_lognormal_balanced.mat');
%r2=load('power_test_lognormal_balanced1.mat');

figure(1); clf;
plot(r1.xs,r1.bs,'b-o'); hold on
plot(r1.xs,r1.bs2,'r-o');  
legend({'2*min(p,1-p)','1-sided abs(score)'});
%plot(r2.xs,r2.bs,'b-o'); hold on
%plot(r2.xs,r2.bs2,'r-o');  

fs=18;
fs1=14;
title('log normal, balanced, snr statistic, n0=20, n1=20','FontSize',fs);
xlabel('\Delta','FontSize',fs);
ylabel('Power','FontSize',fs);
set(gca,'FontSize',fs1);
print -dpdf -f1 lognormal_power_balanced.pdf


fs=18;
fs1=14;
figure(1); clf;
plot(r1.xs,r1.bs,'b-o'); hold on
plot(xs,bs,'r-o');  
legend({'lognormal','normal'});
title('log normal vs. normal, snr, 2*min(p,1-p),  n0=5, n1=50','FontSize',fs);
xlabel('\Delta','FontSize',fs);
ylabel('Power','FontSize',fs);
set(gca,'FontSize',fs1);
print -dpdf -f1 lognormal_vs_normal_unbalanced.pdf




%%%%%
D=read_mit_res_file('~/genepattern/all_aml/all_aml_train.res');
D=read_mit_cls_file(D,'~/genepattern/all_aml/all_aml_train.cls');
[P,rs,P2,S,gp,fwer,fpr]=marker_selection(D.dat,find(D.supdat(1,:)==1),find(D.supdat(2,:)==1),'snr',1000);

[nms,sc,x1,x2,x3,x4,x5]=textread('all_aml_train.txt','%s%f%f%f%f%f%f');
idx=findstrings_list(strvcat(nms),strvcat(D.gacc));
nms=nms(idx);
sc=sc(idx);
x1=x1(idx);
x2=x2(idx);
x3=x3(idx);
x4=x4(idx);
x5=x5(idx);

% finish marker selection comparison!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
cd ~/projects/miRNA/new/TN
load Dx
D=Dx;
[idx,q,p,s,pi0,F]=get_top_markers(Dx,1,struct('method', ...
                                              'ttest_minvar','minvar',(0.75)^2,'nparts_perm',10,'nparts_fix',10),100000,...
                                     struct('method','bonferroni','thresh',0.05));

[ss,si]=sort(-F.s);
save Dx_res.mat idx q p s pi0 F

% look at tmp.mat file in lsfdir

load Dx
load Dx_res.mat
load /xchip/data/gadgetz/lsfres/tmp.mat
[gp,fwer,fpr]=rankbased_pvalues(abs(S),abs(SR));
fdr=calc_fdr_value(F.p);

f=fopen('res.txt','w');
for i=1:length(si)
  i1=si(i);
  fprintf(f,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n',Dx.gacc{i1},F.s(i1),F.p(i1),fpr(i1),fwer(i1),gp(i1),fdr(i1),min(F.p(i1)*length(F.p),1),F.q(i1));
end
fclose(f);

