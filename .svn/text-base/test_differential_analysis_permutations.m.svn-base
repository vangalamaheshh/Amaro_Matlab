clear
cd ~/projects/permutation_tests

D=read_mit_gct_file('norm01.gct');
D=read_mit_cls_file(D,'norm01.cls');

[PR,SR,OP,OS,rs]=differential_analysis_permutations(D,...
                                                  find(D.supdat(1,:)),...
                                                  find(D.supdat(2,:)),'ttest',500,1);


cd /xchip/data/gadgetz/test_perm
% 2n, N(0,1) N(x,1)
% x=0:0.05:1 (20)
% n=5:5:50 (10)
% balanced, unbalanced

N=5000;
Nxv=5; %5;
maxn=32; %8;
res_all=cell(maxn,1);

if (0)
  meta_dat=zeros(Nxv*N,4*maxn);
  for i=1:Nxv
    x=(i-1)*0.5;
    meta_dat((i-1)*N+(1:N),:)=[random('normal',0,1,N,2*maxn) ...
                        random('normal',x,1,N,2*maxn)];
  end
  save meta_dat.mat meta_dat
else
  load meta_dat
end

nvec=[2 4 8 16 32];
use_lsf=0;
clear class
clear l
tic
if use_lsf
  l=lsf('/xchip/data/gadgetz/lsfres/');
  h=zeros(maxn,1);
end
for n=nvec
  if use_lsf
    [l,h(n)]=bsub(l,{'r'},'test_diff_anal_by_parts',{meta_dat,n,maxn,N,Nxv}); ...   
  else
    res_all{n}=test_diff_anal_by_parts(meta_dat,n,maxn,N,Nxv);  
  end
end
if use_lsf
  [l,res]=wait(l); % wait for all

  for n=nvec
    res_all{n}=res{h(n)}.r;
  end
end

tot_time=toc;
save res_all3.mat res_all tot_time

N=2000;
Nxv=5;
maxn=8;
nbin=20;
load res_all1
pvx=cell(maxn,1);

for j=1:maxn
  for i=1:Nxv
    tmp=histc(res_all{j}{1,i},0:1/(nbin):(1+eps));
    pvx{j}(:,i)=tmp(1:nbin);
  end
end


%%% H0: balanced vs. unbalanced

figure(1); clf;
for i=1:maxn
  subplot(2,4,i);
  b1=histc(res_all{i}{1,1},0:1/(nbin):(1+eps));
  b2=histc(res_all{i}{2,1},0:1/(nbin):(1+eps));
  bar([ b1 b2 ],'grouped');
%  comp_dist({res_all{i}{1,1},res_all{i}{2,1}},0.5/nbin:(1/nbin):1);
  ax=axis;
  axis([ -0.05 1.05 ax(3:4)]);
  title(num2str(2*i));
  legend({'unbalanced','balanced'});
end

figure(1); clf;
for i=1:maxn
  subplot(2,4,i);
  b1=histc(res_all{i}{1,1},0:1/(nbin):(1+eps));
  b2=histc(res_all{i}{2,1},0:1/(nbin):(1+eps));
  b1=b1(1:(end-1));
  b2=b2(1:(end-1));
  bar(0.5/nbin:(1/nbin):1,[ b1 b2 ]); %,'grouped');
%  comp_dist({res_all{i}{1,1},res_all{i}{2,1}},0.5/nbin:(1/nbin):1);
  ax=axis;
  axis([ -0.05 1.05 ax(3:4)]);
  title(num2str(2*i));
  legend({'unbalanced','balanced'});
end


save last_point.mat -V6


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=5000;
Nxv=3;
meta_dat2=zeros(Nxv*N,16);
for i=1:Nxv
  x=(i-2);
  meta_dat2((i-1)*N+(1:N),:)=[random('normal',0,1,N,8) ...
                      random('normal',x,1,N,8)];
end
save meta_dat2.mat meta_dat2

D.dat=meta_dat2(10000+(1:5000),:);
% D.dat=meta_dat2((1:5000),:);
cls0=[1:5 9:11];
cls1=[12:16 6:8];

res=cell(2,1);
nparts=3;
verbose('1')
[P,S]=differential_analysis(D,cls0,cls1,'ttest1side');
verbose('2')
[PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1, ...
                                                  'ttest1side',5000,0,nparts,'/xchip/data/gadgetz/lsfres/');
verbose('3')
[Pf,PRf]=fix_Pvalues_by_parts(P,PR,nparts,'/xchip/data/gadgetz/lsfres/');

res{1}=Pf;
verbose('4')

[PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1, ...
                                                  'ttest1side',5000,1,nparts,'/xchip/data/gadgetz/lsfres/');
verbose('5')
[Pf,PRf]=fix_Pvalues_by_parts(P,PR,nparts,'/xchip/data/gadgetz/lsfres/');
verbose('6')

res{2}=Pf;

save res4.mat res -V6




for i=1:5
    figure(i); clf;
    subplot(2,1,1);
    xs=-3:0.01:5;
    plot(xs,normpdf(xs,0,1),'b-','LineWidth',2); hold on
    plot(xs,normpdf(xs,(i-1)*0.5,1),'r-'); hold on
    title(['\Delta = ' num2str((i-1)*0.5)],'FontSize',14);
    
    subplot(2,1,2);
    hist(res_all{4}{1,i},0.5/20:1/20:(1+eps));
    title('histogram of 5000 instances','FontSize',14);
    xlabel('P-value');
end

x=zeros([size(pvx{2}) 5]);
nx=[2 4 8 16 32];
for j=1:5
    x(:,:,j)=pvx{nx(j)};
end
    
for j=[2 4 8 16 32]
    figure(1); clf
    for i=1:5
        sh=subplot(3,2,i);
        hist(res_all{j}{1,i},0.5/20:1/20:(1+eps));
        title(['histogram for \Delta=' num2str((i-1)*0.5) ' (N_s=' num2str(2*j) ')'],'FontSize',10);
        xlabel('P-value');
        shp=get(sh,'Position');
        sh2=axes('position',[shp(1)+0.7*shp(3) shp(2)+0.7*shp(4) shp(3)*0.25 shp(4)*0.25]);
        plot(1:10,1:10);
        box on
    end
    pause
end


figure(1);
for i=1:5
    sh=subplot(2,3,i);
    bar3(squeeze(x(:,i,:)));
    ylabel('P-value','FontSize',14);
    ax=axis;
    axis([ 0.5 5.5 0.5 20.5 0 ax(6)]);
    set(gca,'Ytick',[0.5 5.5 10.5 15.5 20.5]);
    set(gca,'YtickLabel',{'0','0.25','0.5','0.75','1'});
    set(gca,'XTick',1:5);
    set(gca,'XTickLabel',{'4','8','16','32','64'});    
    xlabel('N_s','FontSize',14);
    title(['histograms for \Delta=' num2str((i-1)*0.5) ],'FontSize',14);
end


close all

figure(1); clf;
y=[2 4 8 16 32];
for i=1:5
    subplot(2,3,i);
    comp_dist({res_all{y(i)}{1,1},res_all{y(i)}{2,1}},0.025:0.05:1);
    legend({'unbalanced','balanced'});
    ax=axis;
    axis([-0.025 1.025 0 max(ax(4),350)]);
    title(['N_s=' num2str(2*y(i))],'FontSize',14);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% skewed distrib. equal and unequal class size

cd /xchip/data/gadgetz/aml_all

D=read_mit_res_file('ALL_vs_AML_U95_test.res');
D=read_mit_cls_file(D,'ALL_vs_AML_U95_test.cls');

save ALL_vs_AML_U95_test.mat D

load ALL_vs_AML_U95_test.mat
Dfull=D;

D0=reorder_D_cols(Dfull,find(Dfull.supdat(1,:)));
D1=reorder_D_cols(Dfull,find(Dfull.supdat(2,:)));

D=D0;
D=threshold(D,1);

set1=(D.affy_call==0 | D.tdat==1);
set1a=(D.affy_call==0 & D.tdat>1);
set2=(D.affy_call==2 & D.tdat>1);

figure(1); 
subplot(1,2,1);
hist(log2(D.tdat(set1)),0:0.05:17)
subplot(1,2,1);
hist(log2(D.tdat(set1a)),0:0.05:17)
subplot(1,2,2);
hist(log2(D.tdat(set2)),0:0.05:17)


D=D1;
D=threshold(D,1);

set1=(D.affy_call==0 | D.tdat==1);
set1a=(D.affy_call==0 & D.tdat>1);
set2=(D.affy_call==2 & D.tdat>1);

figure(1); 
subplot(1,2,1);
hist(log2(D.tdat(set1)),0:0.05:17)
subplot(1,2,1);
hist(log2(D.tdat(set1a)),0:0.05:17)
subplot(1,2,2);
hist(log2(D.tdat(set2)),0:0.05:17)





N=1000;
n0=10;
n1=10;


% dat=

