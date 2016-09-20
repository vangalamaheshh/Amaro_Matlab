clear all
close all
cd ~/projects/vogelstein

%% Read CCDS
d=read_dlm_file('/xchip/data/gadgetz/papers/oncomap/supp_data/ccds.txt');
d=cat(1,d{:});
[u,ui,uj]=unique(strvcat(d(:,2)),'rows');

g=d(ui,2);
for i=1:length(ui)
  pos=find(uj==i);
  gsz(i)=sum(str2num(strvcat(d(pos,3))));
end

%% Calc n
% A, C, G, T
% CpG|GpC TpC|GpA A C G T N
ps=[1/4+1/4-1/16-1/32, 1/4-1/32, 6/16]
probs=[ ps(1)*1/4*2 1/4*ps(2)*2 1/4 1/4*ps(3) 1/4*ps(3) 1/4];
n=round(gsz'*probs);
n=[n gsz'];

%% background mutation rates
f_colon=[ 7.73 0.96 0.56 0.95 0.85 0.51 0.55]*1e-6; % colorectal
f_breast=[ 2.99 2.48 0.76 1.38 1.07 0.30 0.55]*1e-6; % breast

f_colon=[ 7.727203 0.9609678 0.5581601 0.9474716 0.8524273 0.5061535 0.552]*1e-6; % colorectal
f_breast=[ 2.985510 2.482500 0.7611274 1.378140 1.065534 0.3036921 0.552]*1e-6; % breast

%% Read CAN genes
clear v
x=read_dlm_file('/xchip/data/gadgetz/papers/oncomap/supp_data/Table S5 - Breast CAN genes 8-29-6_GG.txt');
for i=4:125
  v.p(i-3)=str2num(x{i}{18});
  v.camp(i-3)=str2num(x{i}{5});
  v.symb{i-3}=x{i}{1};
  v.ccds{i-3}=x{i}{2};
end
v_breast=v;

clear v
x=read_dlm_file('/xchip/data/gadgetz/papers/oncomap/supp_data/Table S6 - Colon CAN genes 8-29-6_GG.txt');
for i=3:71
  v.p(i-2)=str2num(x{i}{18});
  v.camp(i-2)=str2num(x{i}{5});
  v.symb{i-2}=x{i}{1};
  v.ccds{i-2}=x{i}{2};
end
v_colon=v;

save before_start.mat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=read_dlm_file(['/xchip/data/gadgetz/papers/oncomap/supp_data/Table S4 -  Somatic mutations identified in breast ' ...
                 'and colorectal cancers 8-29-06_GG.txt']);
for i=5:1678
  m(i-4).ccds=d{i}{1};
  m(i-4).gene=d{i}{2};
  m(i-4).samples=d{i}{4};
  m(i-4).type=d{i}{3};
  m(i-4).mut_genomic=d{i}{5};
  m(i-4).mut_cDNA=d{i}{6};
  m(i-4).mut_aa=d{i}{7};  
end

d=read_dlm_file(['/xchip/data/gadgetz/papers/oncomap/supp_data/Table S3 - Distribution of mutations in individual ' ...
                 'cancers 8-27-06.txt']);

d=read_dlm_file(['/xchip/data/gadgetz/papers/oncomap/supp_data/colon_samples.txt']);                 
for i=1:length(d)
  col_samples(i).name=regexprep(deblank(d{i}{1}),'^[ ]+','');
  col_samples(i).val=regexprep(deblank(d{i}{8}),'^[ ]+','');
end
col_samples(25).name='Hx5';

d=read_dlm_file(['/xchip/data/gadgetz/papers/oncomap/supp_data/breast_samples.txt']);                 
for i=1:length(d)
  breast_samples(i).name=regexprep(deblank(d{i}{1}),'^[ ]+','');
  breast_samples(i).val=regexprep(deblank(d{i}{10}),'^[ ]+','');
end

types={'Breast','Colorectal'};

for t=1:2
  j=strmatch(types{t},{m.type});
  us=cellstr(unique(strvcat({m(j).samples}),'rows'));
  X.sdesc=us;
  X.gacc=g;
  [Mt,m1,m2]=match_string_sets({m.samples},X.sdesc);
  mt=m(m1);
  [Nt,n1,n2]=match_string_sets({mt.gene},X.gacc);
  m2=m2(n1);
  mt2=mt(n1);
  X.ref=cell(length(X.gacc),length(X.sdesc));
  for i=1:length(n1)
    X.ref{n2(i),m2(i)}=[ X.ref{n2(i),m2(i)} mt2(i)];
  end
  X.dat=cellfun('length',X.ref);
  M{t}=X;
end

for t=1:2
  M{1}.sdesc=regexprep(M{1}.sdesc,'[ ]+','');
end

[Mt,m1,m2]=match_string_sets(M{1}.sdesc,{breast_samples.name});
empty_bs=breast_samples(setdiff(1:length(breast_samples),m2));
catn({empty_bs.name});
%1       BB2T 
%2       BB36T
%3       BB37T
%4       BB39T
%5       BB42T
%6       BB44T
M{1}=reorder_D_cols(M{1},m1);
M{1}.supacc='VAL';
M{1}.supdesc='VAL';
v=zeros(1,size(M{1}.dat,2));
v(strmatch('Validation',{breast_samples(m2).val}))=1;
M{1}.supdat=v;

[Mt,m1,m2]=match_string_sets(M{2}.sdesc,{col_samples.name});
empty_cs=col_samples(setdiff(1:length(col_samples),m2));
catn({empty_cs.name});
%1       Mx45
M{2}=reorder_D_cols(M{2},m1);
M{2}.supacc='VAL';
M{2}.supdesc='VAL';
v=zeros(1,size(M{2}.dat,2));
v(strmatch('Validation',{col_samples(m2).val}))=1;
M{2}.supdat=v;

save M.mat M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

types={'Breast','Colon'};
V={};
for i=1:2
  d=read_dlm_file([ '/xchip/data/gadgetz/papers/oncomap/supp_data/' types{i} '_CaMP_GG.txt']);
  c={};
  for j=5:length(d)
    c{j}=d{j}(1:35);
  end
  c=cat(1,c{:});
  tmp=find(~cellfun('isempty',c(:,1)));
  C{i}=c(tmp,:);
  V{i}.sdesc={'N(CG in CpG)','N(G in GpA)','N(C in TpC)','N(A)','N(C)','N(G)','N(T)','N(tot)',...
              'x(CG in CpG)','x(G in GpA|C in TpC)','x(A)','x(C)','x(G)','x(T)','x(IDD)',...
              'mut tot','Discovery bp','Validation bp',...
              'b(1)','b(2)','b(3)','b(4)','b(5)','b(6)','b(7)','b prod'};
  V{i}.gacc=C{i}(:,1);
  V{i}.dat=zeros(size(C{i},1),length([2:9 17:26 28:35]));
  k=0;
  for j=[2:9 17:26 28:35]
    k=k+1;
    V{i}.dat(:,k)=str2num(strvcat(C{i}(:,j)));
  end
end
save V.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determining the cutoff

for i=1:100
  p=unifrnd(0,1,size(n,1),1000);
  sp=sort(p,1);
  q=sp(1,:)*size(p,1);
  t(i)=prctile(q(1,:),10);
  rq=calc_fdr_value(p);
  srq=sort(rq(:));
  t2(i)=srq(100);
  i
end

save ts2.mat t t2
% mean, std = 0.1042, 0.0117

% ==> it is good enough to take a 1000 and set threshold according to 10% percentile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
addpath ~/projects/vogelstein/
cd ~/projects/vogelstein

load before_start.mat

K=35;
[p,vp]=do_analysis('/xchip/data01/gadgetz/vogelstein/colon',n,f_colon,v_colon,K);
svp=sort(vp);
prctile(-log10(svp(1,:)*size(n,1)),90) % 2.391
 

K=35;
[p,vp]=do_analysis('/xchip/data01/gadgetz/vogelstein/breast',n,f_breast,v_breast,K);
svp=sort(vp);
prctile(-log10(svp(1,:)*size(n,1)),90) % 2.3539
 
factor=1.0:0.1:2;
T=1000; %1000
nparts=10; % 10
for i=1:length(factor)
  [p,vp,x,tab]=run_sim(n,K,factor(i)*f_colon,f_colon,T,nparts,1);
  save([ 'colon_' num2str(factor(i)) '_data.mat'],'p','vp'); 
  svp=sort(vp);
  cutoff(i)=prctile(-log10(svp(1,:)*size(n,1)),90);
end
save colon_factor.mat cutoff factor T nparts
figure(1); clf;
plot(factor,cutoff,'o-');
xlabel('Factor');
ylabel('CaMP cutoff');
print_D('colon_factor',{{'epsc'}});

factor=1.0:0.1:2;
T=1000; %1000
nparts=10; % 10
for i=1:length(factor)
  fname=[ 'breast_' num2str(factor(i)) '_data.mat'];
  if exist(fname,'file')
    load(fname);
  else
    [p,vp,x,tab]=run_sim(n,K,factor(i)*f_breast,f_breast,T,nparts,1);
    save([ 'breast_' num2str(factor(i)) '_data.mat'],'p','vp'); 
  end
  svp=sort(vp);
  cutoff(i)=prctile(-log10(svp(1,:)*size(n,1)),90);
end
save breast_factor.mat cutoff factor T nparts
figure(1); clf;
plot(factor,cutoff,'o-');
xlabel('Factor');
ylabel('CaMP cutoff');
print_D('breast_factor',{{'epsc'}});

factor=1.0:0.1:2;
T=1000; %1000
nparts=10; % 10
for i=1:length(factor)
  gf=lognrnd(0,log(factor(i)),size(n,1),1);
  [p,vp,x,tab]=run_sim(n,K,gf*f_colon,f_colon,T,nparts,1);
  save([ 'colon_' num2str(factor(i)) '_gene_var_data.mat'],'p','vp'); 
  svp=sort(vp);
  cutoff(i)=prctile(-log10(svp(1,:)*size(n,1)),90);
end
save colon_gene_var.mat cutoff factor T nparts
figure(1); clf;
plot(log(factor),cutoff,'o-');
xlabel('log std');
ylabel('CaMP cutoff');
print_D('colon_gene_var',{{'epsc'}});

factor=1.0:0.1:2;
T=1000; %1000
nparts=10; % 10
for i=1:length(factor)
  gf=lognrnd(0,log(factor(i)),size(n,1),1);
  [p,vp,x,tab]=run_sim(n,K,gf*f_breast,f_breast,T,nparts,1);
  save([ 'breast_' num2str(factor(i)) '_gene_var_data.mat'],'p','vp'); 
  svp=sort(vp);
  cutoff(i)=prctile(-log10(svp(1,:)*size(n,1)),90);
end
save breast_gene_var.mat cutoff factor T nparts
figure(1); clf;
plot(log(factor),cutoff,'o-');
xlabel('log std');
ylabel('CaMP cutoff');
print_D('breast_gene_var',{{'epsc'}});

%% effect of screening 
run_screening('/xchip/data01/gadgetz/vogelstein/screening/colon',n,f_colon,100,1+(0.05:0.05:1));
run_screening('/xchip/data01/gadgetz/vogelstein/screening/breast',n,f_breast,100,1+(0.05:0.05:1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
types={'Breast','Colorectal'};
load M.mat

M{1}.gsupacc='KNOWN';
M{1}.gsupdesc='KNOWN';
v=zeros(1,size(M{1}.dat,1));
v(grep('^TP53$',M{1}.gacc,1))=1;
M{1}.gsupdat=v;

M{2}.gsupacc='KNOWN';
M{2}.gsupdesc='KNOWN';
v=zeros(1,size(M{2}.dat,1));
v(grep('(^TP53$|^APC$|^KRAS$|^FBXW7$|^SMAD4$)',M{2}.gacc,1))=1;
M{2}.gsupdat=v;

%%%%%%%%%%%%%%%%%%%% DISCOVERY
nnz(any(M{1}.dat(:,find(M{1}.supdat==0)),2)) % 672
nnz(any(M{2}.dat(:,find(M{2}.supdat==0)),2)) % 519
% 1191 genes
sum(sum(M{1}.dat(:,find(M{1}.supdat==0)))) % 732
sum(sum(M{2}.dat(:,find(M{2}.supdat==0)))) % 574
% 1306 mutations

nnz(any(M{1}.dat(~M{1}.gsupdat(1,:),find(M{1}.supdat==0)),2)) % 671
nnz(any(M{2}.dat(~M{2}.gsupdat(1,:),find(M{2}.supdat==0)),2)) % 514

sum(sum(M{1}.dat(~M{1}.gsupdat(1,:),find(M{1}.supdat==0)))) % 722
sum(sum(M{2}.dat(~M{2}.gsupdat(1,:),find(M{2}.supdat==0)))) % 548

%%%%%%%%%%%%%%%%%%%%%%%%%% VALIDATION
nnz(any(M{1}.dat(:,find(M{1}.supdat==1)),2)) % 137
nnz(any(M{2}.dat(:,find(M{2}.supdat==1)),2)) % 105

sum(sum(M{1}.dat(:,find(M{1}.supdat==1)))) % 188
sum(sum(M{2}.dat(:,find(M{2}.supdat==1)))) % 177

nnz(any(M{1}.dat(~M{1}.gsupdat(1,:),find(M{1}.supdat==1)),2)) % 136
nnz(any(M{2}.dat(~M{2}.gsupdat(1,:),find(M{2}.supdat==1)),2)) % 100

sum(sum(M{1}.dat(~M{1}.gsupdat(1,:),find(M{1}.supdat==1)))) % 180
sum(sum(M{2}.dat(~M{2}.gsupdat(1,:),find(M{2}.supdat==1)))) % 130

sum(sum(M{1}.dat(find(M{1}.gsupdat(1,:)),find(M{1}.supdat==1)))) % 8
sum(sum(M{2}.dat(find(M{2}.gsupdat(1,:)),find(M{2}.supdat==1)))) % 47

722/671=1.08; 180/136=1.32
1.32/1.08=1.23

548/514=1.07; 130/100=1.3
1.3/1.07=1.21

x=1+(0.05:0.05:1);
load breast__all_screening
for i=1:length(R)
  r21(i)=R{i}{4}/R{i}{3};
end
[r21(11:12); x(11:12)]

load colon__all_screening
for i=1:length(R)
  r21(i)=R{i}{4}/R{i}{2};
end
[r21(11:12); x(11:12)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=read_dlm_file('/xchip/data/gadgetz/papers/oncomap/CCDS/CCDS.03032005.txt');
d=cat(1,d{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ~/papers/oncomap/supp_data/Table1.mat
c=gene_counts(d,t1);
save c.mat c

tmp=find(~cellfun('isempty',t1(:,1)));
[u,ui,uj]=unique(strvcat(t1(tmp,2)),'rows');
clear same
same(u,M{1}.gacc) % 1
[Mt,m1,m2]=match_string_sets_hash(t1(tmp,2),M{1}.gacc);
ct=zeros(length(ui),size(c,2));
for i=1:length(m1)
  ct(m2(i),:)=ct(m2(i),:)+c(m1(i),:);
end

newn=[ ct(:,1)+ct(:,2) ct(:,3)+ct(:,4) ct(:,5:8) sum(ct,2)];
save newn.mat newn m2

MSAVE=M;
analyze_mutations(1,MSAVE,V,d,f_colon,f_breast,newn);
analyze_mutations(2,MSAVE,V,d,f_colon,f_breast,newn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare sizes in V table to ones I have
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newn(grep('^TP53$',M{1}.gacc,1),:)

k=1;
[Mt,m1,m2]=match_string_sets_hash(M{k}.gacc,V{k}.gacc);
range(m2-(1:length(m2))')

x=V{k}.dat(:,1:8);
x=[ x(:,1) x(:,2)+x(:,3) x(:,4:end)];
xn=newn(m1,:);

y=V{k}.dat(:,17)/length(find(M{k}.supdat==0));
z=V{k}.dat(:,18)/24; % length(find(M{k}.supdat==1));

plot(abs(y-xn(:,end))./xn(:,end))

[sy,si]=sort(y);
figure(1); clf;
plot(y(si),'-'); hold on
plot(z(si),'g-'); hold on
plot(xn(si,end),'r-');
plot(x(si,end),'k-');

median(y./xn(:,end)) 
% 0.9936  --> we can use our n as an approximation for permutations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validation our n by calculating fraction of different bases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range(sum(V{1}.dat(:,1:7),2)-V{1}.dat(:,8)) % make sure sum up correctly

p=V{k}.dat(:,1:7)./repmat(sum(V{1}.dat(:,1:7),2),1,7);
p=[p(:,1) p(:,2)+p(:,3) p(:,4:end)];
pn=xn(:,1:6)./repmat(xn(:,7),1,6);
imagesc(pn-p)
colorbar
zz=[zz(:,1:2) zz(:,3)+zz(:,6) zz(:,4)+zz(5)];
imagesc(zz)
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use our new n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% effect of screening 
addpath ~/projects/vogelstein
run_screening('/xchip/data01/gadgetz/vogelstein/screening/colon',newn,f_colon,100,1+(0:0.05:1),1);
run_screening('/xchip/data01/gadgetz/vogelstein/screening/breast',newn,f_breast,100,1+(0:0.05:1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc corrected p-values for the V data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bf={f_breast,f_colon};
res=calc_corrected_camp(bf,V);
res1=calc_corrected_camp(bf,V,1);
res2poiss=calc_corrected_camp(bf,V,0,2);
res2=calc_corrected_camp(bf,V,0,1);
%save corrected_pv_camp.mat res res1 res2 res2poiss types bf
load MYX
res1_new=calc_corrected_camp(bf,V,1,0,1,MYX);
save corrected_pv_camp_new.mat res res1 res2 res2poiss types bf res1_new
res_005=calc_corrected_camp(bf,V);
save corrected_pv_camp_new_0_001_40.mat res
cbf=bf;
cbf{1}=cbf{1}*2;
cbf{2}=cbf{2}*2;
res1_new_x2=calc_corrected_camp(cbf,V,1,0,1,MYX);


A=res2;
B=res1;
figure(1); clf
subplot(1,2,1);
loglog(A{1}{3},B{1}{3},'.');
axis square
lh=add_x_equ_y_line; set(lh,'color','r');
title('Breast');
xlabel('A'); ylabel('B');
subplot(1,2,2);
loglog(A{2}{3},B{2}{3},'.');
lh=add_x_equ_y_line; set(lh,'color','r');
axis square
title('Colon');
xlabel('A'); ylabel('B');
print_D('compare_p_values_AB',{{'pdf'}});

figure(1); clf
subplot(1,2,1);
loglog(res{1}{3},res1{1}{3},'.');
axis square
lh=add_x_equ_y_line; set(lh,'color','r');
title('Breast');
xlabel('Semi-exact'); ylabel('Exact');
subplot(1,2,2);
loglog(res{2}{3},res1{2}{3},'.');
lh=add_x_equ_y_line; set(lh,'color','r');
axis square
title('Colon');
xlabel('Semi-exact'); ylabel('Exact');
print_D('compare_p_values',{{'pdf'}});

figure(2); clf
subplot(1,2,1);
plot(-log10(res{1}{5}),-log10(res1{1}{5}),'.');
lh=add_x_equ_y_line; set(lh,'color','r');
axis([-2 5 -2 5]);
axis square
xlabel('Semi-exact'); ylabel('Exact');
title('Breast');
subplot(1,2,2);
plot(-log10(res{2}{5}),-log10(res1{2}{5}),'.');
lh=add_x_equ_y_line; set(lh,'color','r');
axis([-2 5 -2 5]);
axis square
title('Colon');
xlabel('Semi-exact'); ylabel('Exact');
print_D('compare_CaMP',{{'pdf'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc corrected p-values for the V data as a function of changing the f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ~/projects/vogelstein
clear all
close all

f_colon=[ 7.727203 0.9609678 0.5581601 0.9474716 0.8524273 0.5061535 0.552]*1e-6; % colorectal
f_breast=[ 2.985510 2.482500 0.7611274 1.378140 1.065534 0.3036921 0.552]*1e-6; % breast
bf={f_breast,f_colon};
save bf_new.mat bf

load V
load corrected_pv_camp.mat
load MYX
load bf_new

[bf{2}(1)  7.7272e-06] %should be 7.7272e-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test s_g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbf=bf;
for k=1:2
  cbf{k}=cbf{k}*1;
end
R=calc_corrected_camp(cbf,V,1,0,0);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 6 28
write_camp_table('s_g_070624',{R},V,1,2);

cbf=bf;
factor=[1.9 1.43];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,0);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 11
write_camp_table('s_g_1.9_1.43_070624',{R},V,{1.9,1.43},2);

cbf=bf;
factor=[2.5 1.9];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,0);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 6

cbf=bf;
factor=[3.31 3.78];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,0);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test s_g^{new}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbf=bf;
for k=1:2
  cbf{k}=cbf{k}*1;
end
R=calc_corrected_camp(cbf,V,1,0,1,MYX);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 34 38
write_camp_table('s_g_NEW_070624',{R},V,1,2);

cbf=bf;
factor=[1.9 1.43];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,1,MYX);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 18
write_camp_table('s_g_NEW_1.9_1.43_070624',{R},V,{1.9,1.43},2);

cbf=bf;
factor=[2.5 1.9];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,1,MYX);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 6

cbf=bf;
factor=[3.31 3.78];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,1,MYX);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test UMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbf=bf;
for k=1:2
  cbf{k}=cbf{k}*1;
end
R=calc_corrected_camp(cbf,V,1,1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 2 11
write_camp_table('UMP_070624',{R},V,1);

cbf=bf;
factor=[1.9 1.43];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 5
write_camp_table('UMP_1.9_1.43_070624',{R},V,1);

cbf=bf;
factor=[2.5 1.9];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 5

cbf=bf;
factor=[3.31 3.78];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test LLRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbf=bf;
for k=1:2
  cbf{k}=cbf{k}*1;
end
R=calc_corrected_camp(cbf,V,1,0,0,[],1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 7 28

cbf=bf;
factor=[1.9 1.43];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,0,[],1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 12

cbf=bf;
factor=[2.5 1.9];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,0,[],1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 8

cbf=bf;
factor=[3.31 3.78];
for k=1:2
  cbf{k}=cbf{k}*factor(k);
end
R=calc_corrected_camp(cbf,V,1,0,0,[],1);
disp([ length(find(R{1}{2}>1)) length(find(R{2}{2}>1))]); % 1 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAT-14
factor=1.0:0.1:2;
nparts=length(factor);
lsfdir='/xchip/data01/gadgetz/vogelstein/lsftmp/';
l=lsf(lsfdir);
h=zeros(nparts,1);
R=cell(length(nparts),1);
for i=1:nparts
  cbf=bf;
  for k=1:2
    cbf{k}=cbf{k}*factor(i);
  end
  R{i}=calc_corrected_camp(cbf,V,1,0,1,MYX); 
end
save corrected_pv_camp_factors_new_cat14_1_0.1_2.mat R factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact (CAT-7)
factor=1.0:0.1:2;
nparts=length(factor);
lsfdir='/xchip/data01/gadgetz/vogelstein/lsftmp/';
l=lsf(lsfdir);
h=zeros(nparts,1);
for i=1:nparts
  cbf=bf;
  for k=1:2
    cbf{k}=cbf{k}*factor(i);
  end
  [l,h(i)]=bsub(l,{'res'},'calc_corrected_camp',{cbf,V,1}); 
end
[l,res]=wait(l); % wait for all
R=cell(length(nparts),1);
for i=1:nparts
  R{i}=res{h(i)}.res;
end  
save corrected_pv_camp_factors_new_exact_1_0.1_2.mat R factor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% semi-exact (CAT-7)
factor=1.0:0.1:2;
%factor=1.0:0.2:5;
nparts=length(factor);
lsfdir='/xchip/data01/gadgetz/vogelstein/lsftmp/';
l=lsf(lsfdir);
h=zeros(nparts,1);
for i=1:nparts
  cbf=bf;
  for k=1:2
    cbf{k}=cbf{k}*factor(i);
  end
  [l,h(i)]=bsub(l,{'res'},'calc_corrected_camp',{cbf,V}); 
%  [l,h(i)]=bsub(l,{'res'},'calc_corrected_camp',{cbf,V}); % for exact add ,1
end
[l,res]=wait(l); % wait for all
R=cell(length(nparts),1);
for i=1:nparts
  R{i}=res{h(i)}.res;
end  
save corrected_pv_camp_factors_new_semi_exact_0_001_25_1_0.1_2.mat R factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% semi-exact (CAT-7) 0.01 25
factor=1.0:0.1:2;
%factor=1.0:0.2:5;
nparts=length(factor);
lsfdir='/xchip/data01/gadgetz/vogelstein/lsftmp/';
l=lsf(lsfdir);
h=zeros(nparts,1);
for i=1:nparts
  cbf=bf;
  for k=1:2
    cbf{k}=cbf{k}*factor(i);
  end
  [l,h(i)]=bsub(l,{'res'},'calc_corrected_camp',{cbf,V,struct('step',0.01,'maxval',25)}); 
%  [l,h(i)]=bsub(l,{'res'},'calc_corrected_camp',{cbf,V}); % for exact add ,1
end
[l,res]=wait(l); % wait for all
R=cell(length(nparts),1);
for i=1:nparts
  R{i}=res{h(i)}.res;
end  
save corrected_pv_camp_factors_new_semi_exact_0_01_25_1_0.1_2.mat R factor
res=calc_corrected_camp(bf,V,struct('step',0.01,'maxval',40));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write tables
factor=1.0:0.1:2;

load corrected_pv_camp_factors_new_cat14_1_0.1_2.mat 
write_camp_table('corrected_camp_factor_new_cat14_1_0.1_2',R,V,factor)

load corrected_pv_camp_factors_new_exact_1_0.1_2.mat
write_camp_table('corrected_camp_factor_new_exact_1_0.1_2',R,V,factor)

load corrected_pv_camp_factors_new_semi_exact_0_001_25_1_0.1_2.mat
write_camp_table('corrected_camp_factor_new_semi_exact_0_001_25_1_0.1_2',R,V,factor)

load corrected_pv_camp_factors_new_semi_exact_0_01_25_1_0.1_2.mat
write_camp_table('corrected_camp_factor_new_semi_exact_0_01_25_1_0.1_2',R,V,factor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_screening('/xchip/data01/gadgetz/vogelstein/screening/colon',newn,f_colon,400,1+(0:0.02:1),1);
run_screening('/xchip/data01/gadgetz/vogelstein/screening/breast',newn,f_breast,400,1+(0:0.02:1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run on many sigmas 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load fact.mat
run_screening('/xchip/data01/gadgetz/vogelstein/screening/breast_s_1.9',newn,f_breast*fact(1),500,1+(0:0.02:1),1);
run_screening('/xchip/data01/gadgetz/vogelstein/screening/colon_s_1.46',newn,f_colon*fact(2),500,1+(0:0.02:1),1);

% run effect 2D 
res=screening_effect_2D('/xchip/data01/gadgetz/vogelstein/screening/colon_2D',...
                        newn,bf{1},bf{1},1.3,1.3,10,1);
lognstd=1:0.05:2;
slognstd=1:0.05:2;

R=run_screening_2D('/xchip/data01/gadgetz/vogelstein/screening/colon_2D',...
                   newn,bf{1},100,lognstd,slognstd,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test effect of sample variability --> NONE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=[];
for jj=1:100
  use_mean=1;
  K=34;
  slognstd=2;
  sfactor=lognrnd(0,log(slognstd),1,K);
  if use_mean
    sfactor=sfactor./mean(sfactor);
  end
  sfactor=[repmat(0.2,1,17) repmat(1.8,1,17)];
  g=1;
  r=0:20;
  for t=1:7
    y=zeros(length(r),K);
    for i=1:K
      y(:,i)=binopdf(r,newn(g,t),bf{1}(t)*sfactor(i))';
    end
    dd=y(:,i);
    for i=2:K
      dd=conv(dd,y(:,i));
    end
    dda=binopdf(r,newn(g,t)*K,bf{1}(t))';
    R(:,t,jj)=dda(1:10)./dd(1:10);
  end
  jj
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

load('/xchip/data01/gadgetz/vogelstein/screening/colon_all_screening.mat');
q=1:0.02:2;
clear X X1
X={}; X1={};
ext='s_1.9_';
for j=1:length(q)
  load([ '/xchip/data01/gadgetz/vogelstein/screening/breast_' ext num2str(q(j)) '_screening.mat']);
  X{1}{j}=sparse(squeeze(sum(x,2)));
  X1{1}{j}=sparse(squeeze(sum(x1,2)));
  j
end
ext='s_1.46_';
for j=1:length(q)
  load([ '/xchip/data01/gadgetz/vogelstein/screening/colon_' ext num2str(q(j)) '_screening.mat']);
  X{2}{j}=sparse(squeeze(sum(x,2)));
  X1{2}{j}=sparse(squeeze(sum(x1,2)));
  j
end
%save XX1.mat X X1
save sXX1.mat X X1

for k=1:2
  nm=[];
  for j=1:length(q)
    nm=[ nm mean(sum(X{k}{j},1))];
  end
  
  nmg=[];
  for j=1:length(q)
    nmg=[ nmg mean(sum(X{k}{j}>0,1))];
  end
  disp(full([ mean(nm) mean(nmg) mean(nm)/mean(nmg) ]));
  
  nm=[];
  for j=1:length(q)
    nm=[ nm mean(sum(X1{k}{j},1))];
  end
  
  nmg=[];
  for j=1:length(q)
    nmg=[ nmg mean(sum(X1{k}{j}>0,1))];
  end
  disp(full([ mean(nm) mean(nmg) mean(nm)/mean(nmg) ]));
end
%   380.0990  369.1357    1.0297
%   46.6114   41.1098    1.1338
%   376.2457  365.6313    1.0290
%   46.0285   40.7260    1.1302
fact(1)=722/380;
fact(2)=548/376;

% simulation II
res=screening_effect('/xchip/data01/gadgetz/vogelstein/screening/breast_f_1.9',newn,fact(1)*bf{1},fact(1)*bf{1}, ...
                     1,500,1);
load breast_f_1.9_1_screening.mat
xx=squeeze(sum(x,2));
disp([ mean(sum(xx>0)) mean(sum(xx)) mean(sum(xx))/mean(sum(xx>0))]);
xx1=squeeze(sum(x1,2));
disp([ mean(sum(xx1>0)) mean(sum(xx1)) mean(sum(xx1))/mean(sum(xx1>0))]);


res=screening_effect('/xchip/data01/gadgetz/vogelstein/screening/colon_f_1.46',newn,fact(2)*bf{2},fact(2)*bf{2},1,500,1);
load colon_f_1.46_1_screening.mat
xx=squeeze(sum(x,2));
disp([ mean(sum(xx>0)) mean(sum(xx)) mean(sum(xx))/mean(sum(xx>0))]);
xx1=squeeze(sum(x1,2));
disp([ mean(sum(xx1>0)) mean(sum(xx1)) mean(sum(xx1))/mean(sum(xx1>0))]);



factor=R{j}{5};
f1=R{j}{6};
f2=R{j}{7};
[ft,fi]=min(abs(f1-median(f1)));
figure(1);
hist(factor,100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power

% 1-(1-0.05)^n=0.99
n=log(0.01)/log(0.95);

load('/xchip/data01/gadgetz/vogelstein/screening/colon_s_1.46_Screening.mat');

lognstd=1:0.02:2;
load('/xchip/data01/gadgetz/vogelstein/screening/colon_s_1.46_all_screening.mat');
ystd={};
ymean={};
for i=1:length(lognstd)
  ymean{2}(i)=mean(R{i}{7});  
  ystd{2}(i)=std(R{i}{7},0,2);
end
load('/xchip/data01/gadgetz/vogelstein/screening/breast_s_1.9_all_screening.mat');
for i=1:length(lognstd)
  ymean{1}(i)=mean(R{i}{7});  
  ystd{1}(i)=std(R{i}{7},0,2);
end

figure(1); clf;
errorbar(log(1:0.02:2),ymean{1}*fact(1),ystd{1}*fact(1),ystd{1}*fact(1)); hold on
errorbar(log(1:0.02:2),ymean{2}*fact(2),ystd{2}*fact(2),ystd{2}*fact(2),'r-'); hold on
fs1=10;
xlabel('std. of log-normal distribution of mutation rates','FontSize',fs1);
ylabel('factor\pmstd','FontSize',fs1);
legend({'Breast','Colon'},'Location','NorthWest');
ax=axis;
axis([-0.01 ax(2) 0.95 ax(4)]);
print_D('mutation_func_of_sig',{{'pdf'},{'fig'}});

figure(2); clf;
col='br';
for k=1:2
  y=cellfun_any('mean(sum(x))',X{k}); 
  y=cat(1,y{:});
  y1=cellfun_any('mean(sum(x))',X1{k}); 
  y1=cat(1,y1{:});
  s=cellfun_any('std(sum(x))',X{k}); 
  s=cat(1,s{:});
  s1=cellfun_any('std(sum(x))',X1{k}); 
  s1=cat(1,s1{:});
  errorbar(log(lognstd),y1,s1,col(k)); hold on
end
ax1=axis;
axis([ -0.01 ax(2) ax1(3:4)]);
ax1=axis;
line(ax(1:2),[180 180],'Color','b');
line(ax(1:2),[130 130],'Color','r');
print_D('screening_mutations_func_sig',{{'pdf'},{'fig'}});

figure(3); clf;
col='br';
for k=1:2
  sz=cellfun_any('sum(x>0)',X{k}); 
  sz=cat(1,sz{:});
  sz1=cellfun_any('sum(x>0)',X1{k}); 
  sz1=cat(1,sz1{:});
  nm=cellfun_any('sum(x)',X{k}); 
  nm=cat(1,nm{:});
  nm1=cellfun_any('sum(x)',X1{k}); 
  nm1=cat(1,nm1{:});
  errorbar(log(lognstd),mean(nm1./sz1,2),std(nm1./sz1,0,2),col(k)); hold on
end
ax1=axis;
axis([ -0.01 ax(2) ax1(3:4)]);
ax1=axis;
line(ax(1:2),[1.32 1.32],'Color','b');
line(ax(1:2),[1.3 1.3],'Color','r');
print_D('mutations_per_genes_func_sig',{{'pdf'},{'fig'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load fact.mat
load newn
load bf

res=screening_effect('/xchip/data01/gadgetz/vogelstein/screening/breast_f_1.9_sig_0.6',newn,fact(1)*bf{1},fact(1)*bf{1}, ...
                     1.8,1,500,1);
save /xchip/data01/gadgetz/vogelstein/screening/breast_f_1.9_sig_0.6_1.8_all.mat res 
cd /xchip/data01/gadgetz/vogelstein/screening/screening

load breast_f_1.9_sig_0.6_1.8_screening.mat
xx=squeeze(sum(x,2));
disp([ mean(sum(xx>0)) mean(sum(xx)) mean(sum(xx)./sum(xx>0))]);
disp([ std(sum(xx>0)) std(sum(xx)) std(sum(xx)./sum(xx>0))]);
xx1=squeeze(sum(x1,2));
disp([ mean(sum(xx1>0)) mean(sum(xx1)) mean(sum(xx1)./sum(xx1>0))]);
disp([ std(sum(xx1>0)) std(sum(xx1)) std(sum(xx1)./sum(xx1>0))]);
disp(  mean(sum(xx1)./sum(xx1>0))/(mean(sum(xx)./sum(xx>0))) );
disp( mean(f2)*fact(2))

% 2.1 too low
res=screening_effect('/xchip/data01/gadgetz/vogelstein/screening/colon_f_1.46_sig_0.74',newn,fact(2)*bf{2},fact(2)*bf{2},2.1,1,500,1);
load colon_f_1.46_sig_0.74_2.1_screening.mat
xx=squeeze(sum(x,2));
disp([ mean(sum(xx>0)) mean(sum(xx)) mean(sum(xx)./sum(xx>0))]);
disp([ std(sum(xx>0)) std(sum(xx)) std(sum(xx)./sum(xx>0))]);
xx1=squeeze(sum(x1,2));
disp([ mean(sum(xx1>0)) mean(sum(xx1)) mean(sum(xx1)./sum(xx1>0))]);
disp([ std(sum(xx1>0)) std(sum(xx1)) std(sum(xx1)./sum(xx1>0))]);
disp(  mean(sum(xx1)./sum(xx1>0))/(mean(sum(xx)./sum(xx>0))) );

% 2.2 good!
res=screening_effect('/xchip/data01/gadgetz/vogelstein/screening/colon_f_1.46_sig_0.79',newn,fact(2)*bf{2},fact(2)*bf{2},2.2,1,500,1);
load colon_f_1.46_sig_0.79_2.2_screening.mat
xx=squeeze(sum(x,2));
disp([ mean(sum(xx>0)) mean(sum(xx)) mean(sum(xx)./sum(xx>0))]);
disp([ std(sum(xx>0)) std(sum(xx)) std(sum(xx)./sum(xx>0))]);
xx1=squeeze(sum(x1,2));
disp([ mean(sum(xx1>0)) mean(sum(xx1)) mean(sum(xx1)./sum(xx1>0))]);
disp([ std(sum(xx1>0)) std(sum(xx1)) std(sum(xx1)./sum(xx1>0))]);
disp(  mean(sum(xx1)./sum(xx1>0))/(mean(sum(xx)./sum(xx>0))) );
disp( mean(f2)*fact(2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% screeing_effect_F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
cd ~/projects/vogelstein
load screening/fact.mat
load newn
load bf


tic
  RES=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_f_1.46_sig_0.79',newn,fact(2)*bf{2},fact(2)*bf{2},2.2,100,1,1,11,24,1);
toc
% 22.3sec /100 w/ tab
% 38.5sec /100 w/o tab
% 3.6sec /20 w/ tab
% 5.8sec /20 w/o tab

% 21.1sec 20 it of 1 w/ tab
% 5.7sec 20 it of 1 w/o tab

% 341 sec 100 it of 20 w/ tab
% 204 sec 1000 it of 1 w/o tab

run_screening_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_s_1.9',newn,bf{1}*fact(1),1000,1,1+(0:0.1:1.5),1,11,24,1);
run_screening_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_s_1.46',newn,bf{2}*fact(2),1000,1,1+(0:0.1:1.5),1,11,24,1);

lognstd=1:0.1:2.5;
load('/xchip/data01/gadgetz/vogelstein/screening_F/breast_s_1.9_all_screening.mat');
RR{1}=R;
load('/xchip/data01/gadgetz/vogelstein/screening_F/colon_s_1.46_all_screening.mat');
RR{2}=R;


ystd={};
ymean={};
for k=1:2
  tmp=cellfun_any('mean(cat(1,x{:,8}))',RR{k});  
  ymean{k}=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,8}))',RR{k});
  ystd{k}=cat(1,tmp{:});
end

figure(1); clf;
errorbar(log(lognstd),ymean{1}*fact(1),ystd{1}*fact(1),ystd{1}*fact(1)); hold on
errorbar(log(lognstd),ymean{2}*fact(2),ystd{2}*fact(2),ystd{2}*fact(2),'r-'); hold on
fs1=10;
xlabel('std. of log-normal distribution of mutation rates','FontSize',fs1);
ylabel('mean scaling factor','FontSize',fs1);
legend({'Breast','Colon'},'Location','NorthWest');
ax=axis;
axis([-0.01 ax(2) 0.95 ax(4)]);
print_D('mutation_F_func_of_sig',{{'pdf'},{'fig'}});

figure(2); clf;
col='br';
for k=1:2
  tmp=cellfun_any('mean(cat(1,x{:,3}))',RR{k});
  y1=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,3}))',RR{k});
  s1=cat(1,tmp{:});
  errorbar(log(lognstd),y1,s1,col(k)); hold on
end
legend({'Breast','Colon'},'Location','NorthWest');
ax1=axis;
axis([ -0.01 ax(2) ax1(3:4)]);
ax1=axis;
line(ax(1:2),[180 180],'Color','b');
line(ax(1:2),[130 130],'Color','r');
print_D('screening_mutations_F_func_sig',{{'pdf'},{'fig'}});

figure(3); clf;
col='br';

y1={};
s1={};
for k=1:2
  for r=1:length(RR{k})
    for i=1:1000
      y1{k}(r,i)=mean(newn(RR{k}{r}{i,6},end));
    end
  end
end

for k=1:2 
  tmp=cellfun_any('mean(cat(1,x{:,9}))',RR{k});
  y1=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,9}))',RR{k});
  s1=cat(1,tmp{:});
  tmp=cellfun_any('mean(cat(1,x{:,10}))',RR{k});
  y2=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,10}))',RR{k});
  s2=cat(1,tmp{:});
  %errorbar(log(lognstd),y1,s1,[col(k) '-']); hold on
  errorbar(log(lognstd),y2,s2,[col(k) ':']); hold on
end
legend({'Breast','Colon'},'Location','NorthWest');
ax1=axis;
axis([ -0.01 ax(2) ax1(3:4)]);
ax1=axis;
line(ax(1:2),[180 180],'Color','b');
line(ax(1:2),[130 130],'Color','r');
print_D('sizes_F_func_sig',{{'pdf'},{'fig'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_stats(RR{1}{1});
run_stats(RR{2}{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulations at 0.57~log(1.77), 0.76~log(2.14)
RES1=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_f_1.9_sig_0.57',newn,fact(1)*bf{1},bf{1},1.77,1000,1,1,11,24,1);
RES2=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_f_1.46_sig_0.76',newn,fact(2)*bf{2},bf{2},2.14,1000,1,1,11,24,1);
run_stats(RES1)
run_stats(RES2)

%save SIM3.mat RES1 RES2
save SIM3V.mat RES1 RES2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulations I again
RES1=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_f_1_sig_0',newn,bf{1},bf{1},1,1000,1,1,11,24,1);
RES2=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_f_1_sig_0',newn,bf{2},bf{2},1,1000,1,1,11,24,1);
run_stats(RES1)
run_stats(RES2)

%save SIM1.mat RES1 RES2
save SIM1V.mat RES1 RES2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulations II again
RES1=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_f_1.9_sig_0',newn,fact(1)*bf{1},bf{1},1,1000,1,1,11,24,1);
RES2=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_f_1.46_sig_0',newn,fact(2)*bf{2},bf{2},1,1000,1,1,11,24,1);
run_stats(RES1)
run_stats(RES2)

save SIM2V.mat RES1 RES2


% calc in screening_effect_F number of CAN genes

% run with the estimated mutation rates 

load V
load M

load est_1
[Mt,m1,m2]=match_string_sets_hash(M{1}.gacc,V{1}.gacc);
dT1=discovery_T1(m1,:);
vT1=validation_T1(m1,:);
my_x=vT1+dT1;
x=V{1}.dat(:,9:15);
imagesc(x-my_x);
MYX{1}=[dT1 vT1];

load est_2
[Mt,m1,m2]=match_string_sets_hash(M{2}.gacc,V{2}.gacc);
dT1=discovery_T1(m1,:);
vT1=validation_T1(m1,:);
my_x=vT1+dT1;
x=V{2}.dat(:,9:15);
imagesc(x-my_x);
MYX{2}=[dT1 vT1];

save MYX.mat MYX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulations I again with new bf
load bf_new
load newn
RES1=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_new_f_1_sig_0',newn,bf{1},bf{1},1,1000,1,1,11,24,1);
RES2=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_new_f_1_sig_0',newn,bf{2},bf{2},1,1000,1,1,11,24,1);
run_stats(RES1)
run_stats(RES2)

%save SIM1.mat RES1 RES2
save SIM1V_new.mat RES1 RES2


%
pv=rand(13023,1000);
spv=sort(pv);
camp=-log10(spv*13023./repmat((1:13023)',1,1000));
ncan=sum(camp>1);
mean(ncan) % 0.108 should be 0.1
std(ncan) % 0.3499 should be 0.3
prctile(ncan,90) % 0 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulations II again w/ new bf
fact(1)=722/380;
fact(2)=548/384;

RES1=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_new_f_1.9_sig_0',newn,fact(1)*bf{1},bf{1},1,1000,1,1,11,24,1);
RES2=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_new_f_1.43_sig_0',newn,fact(2)*bf{2},bf{2},1,1000,1,1,11,24,1);
run_stats(RES1)
run_stats(RES2)
save SIM2V_new.mat RES1 RES2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% screeing_effect_F with new bf and new fact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
cd ~/projects/vogelstein
load new_fact.mat
load newn
load new_bf

run_screening_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_new_s_1.9',newn,bf{1}*fact(1),1000,1,1+(0:0.1:1.5),1,11,24,1);
run_screening_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_new_s_1.43',newn,bf{2}*fact(2),1000,1,1+(0:0.1:1.5),1,11,24,1);

lognstd=1:0.1:2.5;
load('/xchip/data01/gadgetz/vogelstein/screening_F/breast_new_s_1.9_all_screening.mat');
RR{1}=R;
load('/xchip/data01/gadgetz/vogelstein/screening_F/colon_new_s_1.43_all_screening.mat');
RR{2}=R;

save /xchip/data01/gadgetz/vogelstein/screening_F/new_all_screening.mat RR

ystd={};
ymean={};
for k=1:2
  tmp=cellfun_any('mean(cat(1,x{:,8}))',RR{k});  
  ymean{k}=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,8}))',RR{k});
  ystd{k}=cat(1,tmp{:});
end

figure(1); clf;
errorbar(log(lognstd),ymean{1}*fact(1),ystd{1}*fact(1),ystd{1}*fact(1)); hold on
errorbar(log(lognstd),ymean{2}*fact(2),ystd{2}*fact(2),ystd{2}*fact(2),'r-'); hold on
fs1=10;
xlabel('std. of log-normal distribution of mutation rates','FontSize',fs1);
ylabel('mean scaling factor','FontSize',fs1);
legend({'Breast','Colon'},'Location','NorthWest');
ax=axis;
axis([-0.01 ax(2) 0.95 ax(4)]);
print_D('new_mutation_F_func_of_sig',{{'pdf'},{'fig'}});

figure(2); clf;
col='br';
for k=1:2
  tmp=cellfun_any('mean(cat(1,x{:,3}))',RR{k});
  y1=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,3}))',RR{k});
  s1=cat(1,tmp{:});
  errorbar(log(lognstd),y1,s1,col(k)); hold on
end
legend({'Breast','Colon'},'Location','NorthWest');
ax1=axis;
axis([ -0.01 ax(2) ax1(3:4)]);
ax1=axis;
line(ax(1:2),[180 180],'Color','b');
line(ax(1:2),[130 130],'Color','r');
print_D('new_screening_mutations_F_func_sig',{{'pdf'},{'fig'}});

figure(3); clf;
col='br';

y1={};
s1={};
for k=1:2
  for r=1:length(RR{k})
    for i=1:1000
      y1{k}(r,i)=mean(newn(RR{k}{r}{i,6},end));
    end
  end
end

for k=1:2 
  tmp=cellfun_any('mean(cat(1,x{:,9}))',RR{k});
  y1=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,9}))',RR{k});
  s1=cat(1,tmp{:});
  tmp=cellfun_any('mean(cat(1,x{:,10}))',RR{k});
  y2=cat(1,tmp{:});
  tmp=cellfun_any('std(cat(1,x{:,10}))',RR{k});
  s2=cat(1,tmp{:});
  %errorbar(log(lognstd),y1,s1,[col(k) '-']); hold on
  errorbar(log(lognstd),y2,s2,[col(k) ':']); hold on
end
legend({'Breast','Colon'},'Location','NorthWest');
ax1=axis;
axis([ -0.01 ax(2) ax1(3:4)]);
ax1=axis;
line(ax(1:2),[180 180],'Color','b');
line(ax(1:2),[130 130],'Color','r');
print_D('new_sizes_F_func_sig',{{'pdf'},{'fig'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run simulations at 0.58~log(1.79), 0.78~log(2.18)
RES1=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/breast_f_1.9_sig_0.58',newn,fact(1)*bf{1},bf{1},1.79,1000,1,1,11,24,1);
RES2=screening_effect_F('/xchip/data01/gadgetz/vogelstein/screening_F/colon_f_1.43_sig_0.78',newn,fact(2)*bf{2},bf{2},2.18,1000,1,1,11,24,1);
run_stats(RES1)
run_stats(RES2)

save SIM3V_new.mat RES1 RES2

logninv(0.75,0,0.58)/logninv(0.25,0,0.58)
logninv(0.75,0,0.78)/logninv(0.25,0,0.78)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating background mutation rates from the discovery screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0=1.2e-6*1500*3;
r1=r0+0.1;

n1=11;
n2=24;
Ng=13023;
rnk=6;
po=screen_power(r0,r1,n1,n2,Ng,rnk);

r0=1.2e-6*1500*3;
freq=0.1;
r1=r0+freq;

power_figure(0.05,struct('NV',300,'Vstep',3,'ext','.13023'))
power_figure(0.05,struct('NV',300,'Vstep',3,'ext','.13023B'))
power_figure(0.1,struct('NV',200,'Vstep',2))
power_figure(0.2,struct('NV',100,'Vstep',1))

power_figure(0.05,struct('NV',400,'Vstep',4,'ND',100,'Dstep',2,'Ng',3000,'rnk',2,'r0',1500*1.2e-6,'ext','.3000'));
power_figure(0.05,struct('NV',400,'Vstep',4,'ND',100,'Dstep',2,'Ng',6000,'rnk',2,'r0',1500*1.2e-6,'ext','.6000'));


power_figure(0.02,struct('NV',400,'Vstep',4,'ND',100,'Dstep',2,'Ng',3000,'rnk',2,'r0',1500*1.2e-6,'ext','.3000'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power for Pancreatic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
r0=1.2e-6*1500;
r1=r0+0.05;
G=2500;
ext=['.' num2str(G) '.panB'];
[P,C,K]=power_figure(0.05*0.65,struct('NV',150,'Vstep',1,'ND',100,'Dstep',1,'Ng',G,'rnk',5,'r0',r0,'ext',ext));

for i=1:99
  j=100-i;
  cc(i)=C(i,j);
  pp(i)=P(i,j);
end

[mpp,mppi]=max(pp);
figure(2);
plot(pp,'ro-'); hold on;
plot(cc*30/1e6,'bo-');
ax=axis;
line([mppi mppi],ax(3:4),'Color','k');
title(['Power=' num2str(pp(mppi)) '; Cost=' num2str(cc(mppi)*30/1e6) 'M$']);
print_D(['max_power' ext ],{{'fig'},{'pdf'}});
save(['G' ext '.mat'],'P','C','ext','r0','r1');

screen_numbers('pan.12.24.10000.1.2e-6.0.35.1500.0.5.txt',1.2e-6*1500,[0.05 0.1 0.2],0.35,12,24,10000,5,0.5,1500);
screen_numbers('pan.12.24.10000.2.4e-6.0.35.1500.0.5.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,12,24,10000,5,0.5,1500);
screen_numbers('pan.12.24.10000.2.4e-6.0.35.1500.0.1.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,12,24,10000,5,0.1,1500);
screen_numbers('pan.12.24.13000.2.4e-6.0.35.1500.0.1.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,12,24,13000,5,0.1,1500);
screen_numbers('pan.24.48.2500.2.4e-6.0.35.1500.0.5.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,24,48,2500,5,0.5,1500);
screen_numbers('pan.24.72.2500.2.4e-6.0.35.1500.0.5.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,24,72,2500,5,0.5,1500);
screen_numbers('pan.24.72.2500.2.4e-6.0.35.1500.0.1.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,24,72,2500,5,0.1,1500);
screen_numbers('pan.24.96.2500.2.4e-6.0.35.1500.0.1.txt',2*1.2e-6*1500,[0.05 0.1 0.2],0.35,24,96,2500,5,0.1,1500);
%'TCGA.500.3000.1.2e-6.0.35.3.1500.0.1.txt'

single_screen_numbers('TCGA.500.3000.1.2e-6.0.35.3.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,500,3000,3,0.1,1500)
screen_numbers('TCGA.96.404.3000.1.2e-6.0.35.3.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,96,404,3000,3,0.1,1500);
screen_numbers('TCGA.12.84.10000.1.2e-6.0.35.3.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,12,84,10000,3,0.1,1500);

single_screen_numbers('TCGA.500.3000.2.4e-6.0.35.3.1500.0.1.txt',2.4e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,500,3000,3, ...
                      0.1,1500);
screen_numbers('TCGA.96.404.3000.2.4e-6.0.35.3.1500.0.1.txt',2.4e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,96,404,3000,3,0.1,1500);
screen_numbers('TCGA.12.84.10000.2.4e-6.0.35.3.1500.0.1.txt',2.4e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,12,84,10000,3,0.1,1500);

% gene 3000
single_screen_numbers('TCGA.500.3000.1.2e-6.0.35.3.3000.0.1.txt',1.2e-6*3000,[0.01 0.02 0.05 0.1 0.2],0.35,500,3000,3,0.1,3000)
screen_numbers('TCGA.96.404.3000.1.2e-6.0.35.3.3000.0.1.txt',1.2e-6*3000,[0.01 0.02 0.05 0.1 0.2],0.35,96,404,3000,3,0.1,3000);
screen_numbers('TCGA.12.84.10000.1.2e-6.0.35.3.3000.0.1.txt',1.2e-6*3000,[0.01 0.02 0.05 0.1 0.2],0.35,12,84,10000,3,0.1,3000);

% gene 6000
single_screen_numbers('TCGA.500.3000.1.2e-6.0.35.3.6000.0.1.txt',1.2e-6*6000,[0.01 0.02 0.05 0.1 0.2],0.35,500,3000,3,0.1,6000)
screen_numbers('TCGA.96.404.3000.1.2e-6.0.35.3.6000.0.1.txt',1.2e-6*6000,[0.01 0.02 0.05 0.1 0.2],0.35,96,404,3000,3,0.1,6000);
screen_numbers('TCGA.12.84.10000.1.2e-6.0.35.3.6000.0.1.txt',1.2e-6*6000,[0.01 0.02 0.05 0.1 0.2],0.35,12,84,10000,3,0.1,6000);


% gene length 1500, n genes 602
single_screen_numbers('TCGA.500.602.1.2e-6.0.35.1.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,500,602,1,0.1,1500)
screen_numbers('TCGA.96.404.602.1.2e-6.0.35.1.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,96,404,602,1,0.1,1500);

% gene length 1500, n genes 3000
single_screen_numbers('TCGA.500.3000.1.2e-6.0.35.1.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,500,3000,1,0.1,1500)
screen_numbers('TCGA.96.404.3000.1.2e-6.0.35.1.1500.0.1.txt',1.2e-6*1500,[0.01 0.02 0.05 0.1 0.2],0.35,96,404,3000,1,0.1,1500);
[P,C,K,P1,C1,K1]=power_figure(0.05,struct('NV',500,'Vstep',4,'Vplot',10,'ND',100,'Dstep',1,'Dplot',5,'Ng',3000,'rnk',1,'r0',1500*1.2e-6,...
                                          'L',1500,'norm_frac1',0.125,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.35,...
                                          'pointD',94,'pointV',384,'ext','.1500.FLAG0'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test screen_prob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=1e-6*1500;
n1=10;
n2=30;
x=0:20;
p1=poisspdf(x,n1*r);
p2=poisspdf(x,n2*r);
P=p1'*p2;
P(:,1)=0;
P(1,:)=0;
pp=zeros(1,21);
for t=2:20
  for j=1:(t-1)
    i=t-j;
    [ t i j]
    pp(t+1)=pp(t+1)+P(i+1,j+1);
  end
end
v=p1(1)+p2(1)-p1(1)*p2(1);
pp(1)=v;

pp2=screen_prob(x,n1,n2,r,'poisson');
pp2-pp
pp3=screen_prob(x,n1,n2,r,'binom');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MYX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load MYX
for k=1:2
  f=fopen(['X_14_' types{k} '.txt'],'w');
  for j=1:size(MYX{k},1)
    fprintf(f,'%s\t%s\t',V{k}.gacc{j},types{k});
    fprintf(f,[ repmat('%d\t',1,16) '\n' ],[MYX{k}(j,1:7) sum(MYX{k}(j,1:7),2) MYX{k}(j,8:14) sum(MYX{k}(j,8:14),2)]  );
  end
  fclose(f);
end

cd ~/projects/vogelstein/guy
% fid=fopen('242_sym_typ_dis_val_headed.tdt','r');
[l,fid]=read_dlm_file('242_sym_typ_dis_val_headed.tdt',char(9),1);
G=textscan(fid,['%s%s' repmat('%f',1,16)]);
Gdat=cat(2,G{3:end});
fclose(fid);

tps={'Breast','Colorectal'};
for k=1:2
  idx=grep(tps{k},G{2},1);
  GX{k}=Gdat(idx,[1:7 9:15]);
  [Mt,m1,m2]=match_string_sets(G{1}(idx),V{k}.gacc);
  GX{k}=GX{k}(m1,:);
end

%%%%

load ../newn
load ../bf_new
f_breast=bf{1};
f_colon=bf{2};
MSAVE=M;
analyze_mutations(1,MSAVE,V,d,f_colon,f_breast,newn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures for talk in Toronoto, CA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
single_screen_numbers('Toronto.1.txt',1.5e-6*1500,[0.01 0.02 0.03 0.05 0.1 0.2],0.24,500,20000,1,0.1,1500)
screen_numbers('Toronto.2.txt',1.5e-6*1500,[0.01 0.02 0.03 0.05 0.1 0.2],0.24,100,400,20000,1,0.1,1500)

% 0.02
[P,C,K,P1,C1,K1]=power_figure(0.02,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*1.5e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.24,...
                                          'pointD',100,'pointV',400,'point_power',1,'max_power',1,'ext','.Toronoto.2x'));
res02=struct('P',P,'C',C,'K',K,'P1',P1,'C1',C1,'K1',K1);
[P,C,K,P1,C1,K1]=power_figure(0.02,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*1.5e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.24,...
                                          'pointD',100,'pointV',400,'point_power',0,'max_power',0,'ext','.Toronoto.2x.no_max_power'));

% 0.03
[P,C,K,P1,C1,K1]=power_figure(0.03,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*1.5e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.24,...
                                          'pointD',100,'pointV',400,'point_power',1,'max_power',1,'ext','.Toronoto.2xS'));
res03=struct('P',P,'C',C,'K',K,'P1',P1,'C1',C1,'K1',K1);
[P,C,K,P1,C1,K1]=power_figure(0.03,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*1.5e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.24,...
                                          'pointD',100,'pointV',400,'point_power',0,'max_power',0,'ext','.Toronoto.2xS.no_max_power'));

%0.05
[P,C,K,P1,C1,K1]=power_figure(0.05,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*1.5e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.24,...
                                          'pointD',100,'pointV',400,'point_power',1,'max_power',1,'ext','.Toronoto.2x'));
res05=struct('P',P,'C',C,'K',K,'P1',P1,'C1',C1,'K1',K1);
[P,C,K,P1,C1,K1]=power_figure(0.05,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*1.5e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.24,...
                                          'pointD',100,'pointV',400,'point_power',0,'max_power',0,'ext','.Toronoto.2x.no_max_power'));

qq=zeros(150,1);
for i=1:150
  qq(i)=res05.P(i,floor((500-i)/5)+1);
  plot(floor((500-i)/5)+1,i,'wv');

end

figure(1);
clf;
fnr=0.24;
F=1-fnr;
r0=1.5e-6*1500;
mut_rate=[0.01 0.02 0.03 0.05 0.1 0.2];
r1=r0+(mut_rate*F);
L=1500;
Ng=20000;
rnk=1;
po=zeros(length(r1),500);
ks=zeros(length(r1),500);
q=0.1;
for n=1:500
  for i=1:length(r1)
    [po(i,n),ks(i,n)]=single_screen_power(r0,r1(i),n,Ng,rnk,q,L);
  end
  if mod(n,10)==0
    disp(n);
  end
end

po_sm=zeros(length(r1),500);
q=0.1;
for n=1:500
  for i=1:length(r1)
    [po_sm(i,n),ks(i,n)]=single_screen_power_smooth(r0,r1(i),n,Ng,rnk,q,L);
  end
  if mod(n,10)==0
    disp(n);
  end
end

figure(1); clf;
fs=14;
%plot(po(2:(end-1),:)','LineWidth',3); hold on
plot(po_sm(2:(end-1),:)','LineWidth',3); %,'LineStyle',':'); hold on
legend({'2% (1.5%)','3% (2.3%)','5% (3.8%)','10% (7.6%)'},'Location','northwest');
xlabel('No. of samples','FontSize',fs);
ylabel('Power','FontSize',fs);
set(gca,'FontSize',11);
print_D('power_single_smooth',{{'fig'},{'pdf'},{'png','-r180'}});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare to Parmigiani's numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ~/projects/vogelstein/Parmigiani
load ../MYX
load ../V
load ../M
load ../bf_new
load ../newn

d=read_table('CCDSsizes06.txt','%s%f%f%f%f%f%f%f',char(9),1);
d=cat(1,d{:});
d.x=cat(2,d.dat{2:end});
clear D
D.gacc=d.dat{1};
D.gdesc=D.gacc;
D.sdesc=d.headers{1}(2:end);
D.dat=d.x;

[Mt,m1,m2]=match_string_sets_hash(D.gacc,g);

P_g=g(m2);
P_n=newn(m2,:);
D=reorder_D_rows(D,m1);

%% 294 % CCDS5229.1 % AKAP12
qq=fastaread('CCDS_nucleotide.20070227.fna','BLOCKREAD',4915);
D.dat(294,:)
P_n(294,:)

s=['N' qq.Sequence 'N'];
cc=0;
dd=0;
for i=2:(length(s)-1)
  if s(i)=='C' & s(i-1)~='T' & s(i+1)~='G'
    cc=cc+1;
    if s(i)=='C' & s(i-1)=='G'
      dd=dd+1;
    end
    fprintf(1,'%d\t%d\t%s\n',cc,dd,s(i-1:i+1));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ~/projects/power
single_screen_numbers('Toronto.1.txt',3e-6*1500,[0.01 0.02 0.03 0.05 0.1 0.2],0.2,500,20000,1,0.1,1500)
screen_numbers('Toronto.2.txt',3e-6*1500,[0.01 0.02 0.03 0.05 0.1 0.2],0.2,100,400,20000,1,0.1,1500)

% 0.02
[P,C,K,P1,C1,K1]=power_figure(0.02,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*3e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.2,...
                                          'pointD',100,'pointV',400,'point_power',1,'max_power',1,'ext','.Toronoto.2x'));
res02=struct('P',P,'C',C,'K',K,'P1',P1,'C1',C1,'K1',K1);
[P,C,K,P1,C1,K1]=power_figure(0.02,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*3e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.2,...
                                          'pointD',100,'pointV',400,'point_power',0,'max_power',0,'ext','.Toronoto.2x.no_max_power'));

% 0.03
[P,C,K,P1,C1,K1]=power_figure(0.03,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*3e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.2,...
                                          'pointD',100,'pointV',400,'point_power',1,'max_power',1,'ext','.Toronoto.2xS'));
res03=struct('P',P,'C',C,'K',K,'P1',P1,'C1',C1,'K1',K1);
[P,C,K,P1,C1,K1]=power_figure(0.03,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*3e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.2,...
                                          'pointD',100,'pointV',400,'point_power',0,'max_power',0,'ext','.Toronoto.2xS.no_max_power'));

%0.05
[P,C,K,P1,C1,K1]=power_figure(0.05,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*3e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.2,...
                                          'pointD',100,'pointV',400,'point_power',1,'max_power',1,'ext','.Toronoto.2x'));
res05=struct('P',P,'C',C,'K',K,'P1',P1,'C1',C1,'K1',K1);
[P,C,K,P1,C1,K1]=power_figure(0.05,struct('NV',500,'Vstep',4,'Vplot',10,'ND',150,'Dstep',1,'Dplot',5,'Ng',20000,'rnk',1,'r0',1500*3e-6,...
                                          'L',1500,'norm_frac1',1,'norm_frac2',1,'cost_per_amp',1.75,'fnr',0.2,...
                                          'pointD',100,'pointV',400,'point_power',0,'max_power',0,'ext','.Toronoto.2x.no_max_power'));

qq=zeros(150,1);
for i=1:150
  qq(i)=res05.P(i,floor((500-i)/5)+1);
  plot(floor((500-i)/5)+1,i,'wv');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnr=0.2;
F=1-fnr;
L=1500;
r0=3e-6*L;
%mut_rate=[0.01 0.02 0.03 0.05 0.1 0.2];
mut_rate=0.005:0.005:1;
r1=r0+(mut_rate*F);
Ng=20000;
rnk=1;
Ns=4000;
Ns_step=20;
po=zeros(length(r1),Ns/Ns_step);
ks=zeros(length(r1),Ns/Ns_step);

if (0)
  q=0.1;
  for n=1:Ns
    for i=1:length(r1)
      [po(i,n),ks(i,n)]=single_screen_power(r0,r1(i),n,Ng,rnk,q,L);
    end
    if mod(n,10)==0
      disp(n);
    end
  end
end

po_sm=zeros(length(r1),Ns/Ns_step);
q=0.1;
for n=Ns_step:Ns_step:Ns
  for i=1:length(r1)
    [po_sm(i,round(n/Ns_step)),ks(i,round(n/Ns_step))]=single_screen_power_smooth(r0,r1(i),n,Ng,rnk,q,L);
  end
  if mod(n,10)==0
    disp(n);
  end
end

figure(1); clf;
fs=14;
%plot(po(2:(end-1),:)','LineWidth',3); hold on
idx=[2 4 6 10 20 40];
plot(po_sm(idx(2:(end-1)),:)','LineWidth',3); %,'LineStyle',':'); hold on
legend({'2% (1.5%)','3% (2.3%)','5% (3.8%)','10% (7.6%)'},'Location','northwest');
xlabel('No. of samples','FontSize',fs);
ylabel('Power','FontSize',fs);
set(gca,'FontSize',11);
print_D('power_single_smooth',{{'fig'},{'pdf'},{'png','-r180'}});

n1=100;
n2=400;
x=0:50;
%i=3;
%p=screen_prob(x,n1,n2,r1(i),'binom',L,r0);

pv_cutoff=0.15;
force_ks=[];
Ng=1325;
rnk=1;
q=0.1;
i=3;
[po,ks]=single_screen_power_smooth(r0,r1(i),n1+n2,Ng,rnk,q,L);
[po,ks]=screen_power_smooth(r0,r1(i),n1,n2,Ng,rnk,q,L,force_ks,pv_cutoff);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final figures for paper 2008/11/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1: Cartoon

% Figure 2
% Effect of number of samples on statistical power to detect cancer-associated genes. (a) Distribution of number
% of mutations in a hypothetical cancer gene which is mutated in 3% of patients (purple bars) compared to a gene
% unrelated to cancer and having only background mutations (yellow bars). Calculation for 20, 100 and 500 samples
% assumes a typical gene length of 1500 bases, background mutation rate of 3x10-6 mutations/base and missing data
% rate of 20%. Vertical line represents the minimal number of mutations needed to reach genome-wide significance,
% p-value < 0.1/20,000. The power to detect the cancer gene is the probability of observing at least this many
% mutations in the sequencing experiment. (b) Power to detect cancer-associated genes mutated with varied
% prevalence (2%, 3%, 5% and 10% of patients) as a function of number of samples. Points represent power
% calculated in (a).  (c) Number of samples needed to reach 80% power as a function of mutation prevalence in
% patients.

clear all 
close all

fnr=0.2;
F=1-fnr;
L=1500;
r0=3e-6*L;
mut_rate=0.03;
r1=r0+(mut_rate*F);
pv=0.1/20e3;

figure(1); clf;
gr=make_subplotgrid(2,[1 1 1],[1 1],[1 1 1 1 0.5],0.2,0.35);
Ns=[20 100 500];
max_mut=20;
for i=1:3
  subplotgrid(gr,i,1);
  bg_mut=binopdf(0:max_mut,Ns(i)*L,r0/L);
  cancer_mut=binopdf(0:max_mut,Ns(i)*L,r1/L);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns(i),20000,1,0.1,L)
  bh=bar(0:max_mut,[bg_mut; cancer_mut]','grouped');
  set(bh(1),'FaceColor','cyan');
  set(bh(2),'FaceColor','red');
  ax=axis;
  axis([-0.7 max_mut+0.7 ax(3:4)]);
  box off
  title([ num2str(Ns(i)) ' Samples'],'FontSize',12);
  ylabel('Probability');
  if (i==3)
    xlabel('No. of mutations');
  end
  tail=1-cumsum([0 bg_mut]);
  cancer_tail=1-cumsum([0 cancer_mut]);
  cutoff=min(find(tail<=pv))-1;
  line([ cutoff-0.5 cutoff-0.5],[ax(3) 0.8*ax(4)],'LineStyle','-');
  line([ cutoff-0.5 cutoff+5 cutoff+5-0.5],[0.8*ax(4) 0.8*ax(4) 0.82*ax(4) ])
  text(cutoff,ax(4)*0.75,['Power = ' num2str(cancer_tail(cutoff+2),'%4.2f')]);
  text(cutoff,ax(4)*0.85,['Significance=1/200,000']);
end

print_D('Figure2a',{{'fig'},{'pdf'},{'png','-r180'}});

figure(2); clf;
%powers=[20 50 80 90];
powers=[80];
for i=1:length(powers)
  lh=loglog(100*(0.005:0.005:1),sum(po_sm<powers(i)/100,2)*Ns_step,'k.-'); hold on;
  set(lh,'LineWidth',1.5);
end

z_a=norminv(1-pv,0,1);
z_p=norminv(0.2,0,1);
t=(r1-r0);
n_est1=(z_p)^2./t;
n_est2=(z_p^2+z_a^2)*r0./(t.^2);
n_est3=-2*z_p*z_a*r0*sqrt(1+t./r0)./(t.^2);


n_est=(z_p)^2./t+(z_p^2+z_a^2)*r0./(t.^2)-2*z_p*z_a*r0*sqrt(1+t./r0)./(t.^2);

lh_est=loglog(100*(0.005:0.005:1),n_est,'r-');
lh_est_part(1)=loglog(100*(0.005:0.005:1),n_est1,'r:');
lh_est_part(2)=loglog(100*(0.005:0.005:1),n_est2,'g:');
lh_est_part(3)=loglog(100*(0.005:0.005:1),n_est3,'b:');


xlabel('Mutation prevalence [%]');
ylabel('No. of samples needed');
ylim([1 2100]);
xlim([0.1 100]);
print_D('Figure2c',{{'fig'},{'pdf'},{'png','-r180'}});

[po,ks]=single_screen_power_smooth(r0,0.01,4000,Ng,rnk,0.1,L);


figure(3); clf;
xcoord=0:100;
semilogy(xcoord,1-normcdf(xcoord,18,sqrt(18))); hold on
semilogy(xcoord,1-binocdf(xcoord,4e3*1.5e3,3e-6),'r-'); hold on
semilogy(xcoord,1-poiss_norm_approx_cdf(xcoord,18),'g-'); hold on
ax=axis;
line(ax(1:2),[5e-6 5e-6]);

figure(3); clf;
xcoord=0:100;
semilogy(xcoord,normcdf(xcoord+0.5,18,sqrt(18))-normcdf(xcoord-0.5,18,sqrt(18))); hold on
semilogy(xcoord,binopdf(xcoord,4e3*1.5e3,3e-6),'r-'); hold on
semilogy(xcoord,poiss_norm_approx_pdf(xcoord,18),'g-'); hold on
ax=axis;
line(ax(1:2),[5e-6 5e-6]);

figure(1); clf;
xcoord=100:-1:0;
tmp1=poiss_norm_approx_pdf(xcoord,18);
tmp2=poiss_norm_approx_cdf(xcoord,18);
loglog(cumsum(tmp1),tmp2,'.');



% Figure 3
%Effect of two-stage experimental design on power. Y-axis represents the number of samples in the discovery phase
%and X-axis represents the number of samples in the validation phase. Color represents the power to detect a
%cancer gene mutated in 3% of patients. Black lines connect experiments with similar amount of sequencing
%(iso-sequencing lines) which roughly represent experimental cost. Dashed line connects experiments which yield
%the maximal power for their sequencing capacity.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation for the Melanoa BRAF inhibitor P01 2011/1/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnr=0.23;
F=1-fnr;
L=1500;
r0=2e-6*L;
mut_rate=0.02;
r1=r0+(mut_rate*F);
pv=0.1/20e3;
for i=1:100
  for j=1:100
    Ns=i;
    mut_rate=j/100;
    r1=r0+(mut_rate*F);
    [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
    P(i,j)=po;
  end
  disp(i);
end

imagesc(P>=0.8);

x=[];
for i=1:100
  tmp=find(P(:,i)>=0.8);
  if isempty(tmp)
    x(i)=NaN;
  else
    x(i)=tmp(1);
  end
end

figure(1); clf;
h=plot(x);
grid on
set(h,'LineWidth',2);
xlabel('% of relapses with mutation','FontSize',18);
ylabel('Number of patients needed','FontSize',18);
set(gca,'FontSize',16);
hold on;
plot(14,60,'.','MarkerSize',40,'Color','r');

cd('/xchip/data01/gadgetz/power');
print_D('Melanoma_relapse_power',{{'pdf'}});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation for NHGRI renewal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all
cd('/Users/gadgetz/Google Drive/_DROPBOX/work/matlab/gaddy/power/');


fnr=0.23;
F=1-fnr;
L=1500;
r0=1.5e-6*L;
mut_rate=0.03;
r1=r0+(mut_rate*F);
pv=0.1/20e3;
Ns=400;
% [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
r0s=logspace(log10(0.5e-6),log10(20e-6),90);
xv=zeros(90,1);
prev=log(200);
ug=0;
for i=1:length(xv)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),20000,1,0.1,L,5000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,20000,1,0.1,L,5000,ug))-0.85,prev);
  xv(i)=exp(x);
  prev=x;
end

if (0)
  % use xv as initialization
  xv2=zeros(90,1);
  prev=log(200);
  ug=1;
  for i=1:length(xv2)
    r0=r0s(i)*L;
    r1=r0+(mut_rate*F);
    x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),20000,1,0.1,L,5000,ug)+ ...
                   (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,20000,1,0.1,L,5000,ug))-0.85,log(xv(i)));
    xv2(i)=exp(x);
  end
end

xv1=zeros(90,1);
prev=log(200);
ug=0;
mut_rate=0.02;
for i=1:length(xv1)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),20000,1,0.1,L,5000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,20000,1,0.1,L,5000,ug))-0.85,prev);
  xv1(i)=exp(x);
  prev=x;
end

xv2=zeros(90,1);
prev=log(200);
ug=0;
mut_rate=0.01;
for i=1:length(xv2)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),20000,1,0.1,L,15000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,20000,1,0.1,L,15000,ug))-0.85,prev);
  xv2(i)=exp(x);
  prev=x;
end


xv3=zeros(90,1);
prev=log(200);
ug=0;
mut_rate=0.05;
for i=1:length(xv2)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),20000,1,0.1,L,15000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,20000,1,0.1,L,15000,ug))-0.85,prev);
  xv3(i)=exp(x);
  prev=x;
end

figure(1); clf;
lh(4)=loglog(r0s,xv3,'k');
hold on
lh(1)=loglog(r0s,xv,'b');
lh(2)=loglog(r0s,xv1,'g');
lh(3)=loglog(r0s,xv2,'r');
xlim([0.5e-6 20e-6]);
ylim([100 2e4]);
grid on
set(lh,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'XTick',1e-6*[0.5 1:10 20],'XTickLabel',{'    0.5','1','2','3','4','5','','','','','10', ...
                    '20'});
set(gca,'YTick',[100:100:1000 2000:1000:10000 20000],'YTickLabel',{'100','200','300','400','500','','','','',...
                    '1000','2000','3000','4000','5000','','','','','10,000','20,000'});
legend({'5%','3%','2%','1%'},'Location','NorthWest');
xlabel('Mutation rate [mut/Mb]','FontSize',18);
ylabel('Number of patients for power of 85%','FontSize',18);

% cd ~/projects/power
print_D('Number_of_samples_vs_mut_rate_5_to_1',{{'pdf'},{'png','-r180'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation for Matthew in lung adeno cancer (LUAD rate ~ 12 mut/Mb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
fnr=0.23;
F=1-fnr;
L=1500;

for j=1:12
  disp(j);
  r0=j*1e-6*L; 
  mut_rate=0.03;
  r1=r0+(mut_rate*F);
  pv=0.1/20e3;
  % Ns=400;
  % [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  Nss=10:10:1000;
  xv=zeros(length(Nss),1);
  prev=-log(0.01);
  ug=0;
  opt=optimset('Display','off');
  for i=length(xv):-1:1
    Ns=Nss(i);
    x=fsolve(@(x) single_screen_power_smooth(r0,r0+(exp(-x)*F),Ns,20000,1,0.1,L,5000,ug)-0.85,prev,opt);
    xv(i)=exp(-x);
    prev=x;
  end
  xvs{j}=xv;
end

figure(1); clf;
for j=1:12
  loglog(10:10:1000,xvs{j}'); hold on;
end

figure(1); clf;
%for j=1:12
  loglog(10:10:1000,xvs{12}'); hold on;
%end
grid on
set(gca,'FontSize',14);
xlabel('Number of patients','FontSize',16);
ylabel('Freq. with power \geq 85%','FontSize',16);
title('Mutation rate = 12 mut/Mb');
cd ~/projects/power
print_D('Freq_with_85_power_vs_number_of_patients.12mut_per_mb',{{'pdf'},{'png','-r180'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation for Matthew in lung adeno cancer (LUAD rate ~ 12 mut/Mb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
fnr=0.23;
F=1-fnr;
L=1500;

r0=12e-6*L; 
mut_rates=0.01:0.01:0.3;
pv=0.1/20e3;
Ns=65;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv1=xv;

Ns=272;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv2=xv;

Ns=152;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv152=xv;

Ns=172;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv172=xv;


Ns=152+207;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv359=xv;

Ns=172+207;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv379=xv;

figure(1); clf;
plot(mut_rates*100,xv152,'b'); hold on;
plot(mut_rates*100,xv172,'b--'); hold on;
plot(mut_rates*100,xv359,'r'); hold on;
plot(mut_rates*100,xv379,'r--'); hold on;

legend({'152 patients','172 patients','359 patients','379 patients'},'Location','SouthEast');
xlabel('% of patients with mutations','FontSize',16);
ylabel('Power','FontSize',16);
set(gca,'FontSize',14);
title('BMR = 12 mut/Mb','FontSize',18);
grid on

cd ~/projects/power
print_D('Power_vs_freq_of_mut_for_152_172_and_359_379_patients.12mut_per_mb',{{'pdf'},{'png','-r180'}});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation for Nir in melanoma (MEL rate ~ 12 mut/Mb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fnr=0.23;
F=1-fnr;
L=1500;
r0=12e-6*L;
mut_rates=0.01:0.01:0.3;
% r1=r0+(mut_rate*F);
pv=0.1/20e3;
Nss=10:10:100;
% [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);

xv=zeros(length(mut_rates),length(Nss));
for i=1:length(xv)
  for j=1:length(Nss)
    mut_rate=mut_rates(i);
    r1=r0+(mut_rate*F);
    Ns=Nss(j);
    [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
    xv(i,j)=po;
  end
end

figure(1); clf;
imagesc(xv);
set(gca,'XTick',1:11,'XTickLabel',cellstr(num2str((10:10:100)')));
set(gca,'YTick',1:30,'YTickLabel',cellstr(num2str((1:30)')));
caxis([0 1]);
colorbar
title('Power');
xlabel('No. of samples');
ylabel('Percent patients with mutation');


cd ~/projects/power
print_D('Power.12mut_per_mb',{{'pdf'},{'png','-r180'}});


%% calc power for Levi's grant PLCO
%%%------------------------------------------------------
clear all
close all

fnr=0.23;
F=1-fnr;
L=1500;

r0=10e-6*L %used to be  6e-6*L; % changed to 10 to reply to reviewer
mut_rates=0.01:0.01:0.4;
pv=0.1/20e3;

Ns=275;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv275=xv;

Ns=30;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv30=xv;

Ns=80;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv80=xv;

Ns=385;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv385=xv;

figure(1); clf;
plot(mut_rates*100,xv30,'r'); hold on;
plot(mut_rates*100,xv80,'g'); hold on;
plot(mut_rates*100,xv275,'b'); hold on;
plot(mut_rates*100,xv385,'k'); hold on;

legend({'30 CRC missed patients','80 CRC screening-detected patients','275 CRC usual care patients','385 all together'},'Location','SouthEast');
xlabel('% of patients with mutations','FontSize',16);
ylabel('Power','FontSize',16);
set(gca,'FontSize',14);
title('BMR = 10 mut/Mb','FontSize',18);
grid on


%% Levi's grant power for chi-square as function of effect size (odds ratio)

nx=99;
nnotx=198;
tot=nx+nnotx;

px=nx/tot;
pnotx=1-px;

py=0.2;
ny=tot*py;
nnoty=tot-ny;
pnoty=1-py;

R=2;





odsr=10.^([-1:0.01:1]);

p=zeros(length(odsr),1);

for i=1:length(odsr)
    if abs(odsr(i)-1)<eps
        p(i)=1;
    else
        R=odsr(i);
        s=sqrt((1+(px+py)*(R-1)).^2+4*R*(1-R)*px*py);
        m(2,2)=tot*(1+(px+py)*(R-1)-s)/(2*(R-1));
        m(1,2)=ny-m(2,2);
        m(2,1)=nx-m(2,2);
        m(1,1)=tot-m(1,2)-m(2,1)-m(2,2);
        
        
        e=tot*[pnotx; px]*[pnoty py];
        score=sum(sum((m-e).^2./e,1),2);
        
        disp(odsr(i));
        disp([ m e]);
        disp(score);
        disp('---------');
       
        p(i)=1-chi2cdf(score,1);
    end
end

figure(1); clf;
loglog(odsr,p);
set(gca,'FontSize',12);
ylim([.001,1]);

title([num2str(nx) ' vs. ' num2str(nnotx) '; mutation freq.=' num2str(py*100) '%'],'FontSize',16);
xlabel('Odds Ratio');
ylabel('P-value');

odsr5=odsr;
p5=p;


%%


figure(1); clf;
loglog(odsr5,p5,'r'); hold on;
loglog(odsr10,p10,'g'); hold on;
loglog(odsr15,p15,'b'); hold on;
loglog(odsr20,p20,'k'); hold on;

set(gca,'FontSize',12);
ylim([.001,1]);

title([num2str(nx) ' vs. ' num2str(nnotx) ],'FontSize',16);
xlabel('Odds Ratio');
ylabel('P-value');
ax=axis;
line(ax(1:2),[0.05 0.05]);
legend('5%','10%','15%','20%');
grid on

cd ('/Users/gadgetz/Google Drive/_DROPBOX/work/grants/PLCO_2013');
saveas(gcf,['power_for_odds_ratio.' num2str(nx) ' vs ' num2str(nnotx) '.png'],'png');

%% calc power for Grant with Adam Bass, Andy Chan PLCO
%%%------------------------------------------------------
clear all
close all

fnr=0.23;
F=1-fnr;
L=1500;

r0=10e-6*L %used to be  6e-6*L; % changed to 10 to reply to reviewer
mut_rates=0.01:0.001:0.25;
pv=0.1/20e3;

Ns=298;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv298=xv;

Ns=97;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv97=xv;

Ns=201;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
  xv(i)=po;
end
xv201=xv;

figure(1); clf;
plot(mut_rates*100,xv97,'r'); hold on;
plot(mut_rates*100,xv201,'g'); hold on;
plot(mut_rates*100,xv298,'b'); hold on;

legend({'97 interval CRCs','201 screen-detected CRCs','298 total cases'},'Location','SouthEast');
xlabel('% of patients with mutations','FontSize',16);
ylabel('Power','FontSize',16);
set(gca,'FontSize',14);
title('BMR = 10 mut/Mb','FontSize',18);
grid on

p80i_97=find(xv97>=0.8,1,'first'); disp([ 97 mut_rates(p80i_97)*100]);
p80i_201=find(xv201>=0.8,1,'first'); disp([ 201 mut_rates(p80i_201)*100]);
p80i_298=find(xv298>=0.8,1,'first'); disp([ 298 mut_rates(p80i_298)*100]);

%%
output_dir='/Users/gadgetz/Google Drive/_DROPBOX/work/grants/PLCO_2013/';
saveas(gcf,[output_dir 'power_detctoin_of_cancer_genes.BMR10.97.201.298.png'],'png');



%% calc power for WGS --- DID NOT WORK
%%%------------------------------------------------------
clear all 
close all
cd('/Users/gadgetz/Google Drive/_DROPBOX/work/matlab/gaddy/power/');


fnr=0.23;
F=1-fnr;
L=130;
r0=1.5e-6*L;
mut_rate=0.03;
r1=r0+(mut_rate*F);
pv=0.1/2375714; % 20e3;
NNN=2375714;
Ns=400;
% [po,ks]=single_screen_power_smooth(r0,r1,Ns,20000,1,0.1,L);
r0s=logspace(log10(0.5e-6),log10(20e-6),90);
xv=zeros(90,1);
prev=log(200);
ug=0;
for i=1:length(xv)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),NNN,1,0.1,L,15000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,NNN,1,0.1,L,15000,ug))-0.85,prev);
  xv(i)=exp(x);
  prev=x;
end


xv1=zeros(90,1);
prev=log(200);
ug=0;
mut_rate=0.15;
for i=1:length(xv1)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),NNN,1,0.1,L,15000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,NNN,1,0.1,L,15000,ug))-0.85,prev);
  xv1(i)=exp(x);
  prev=x;
end

xv2=zeros(90,1);
prev=log(200);
ug=0;
mut_rate=0.10;
for i=1:length(xv2)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),NNN,1,0.1,L,15000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,NNN,1,0.1,L,15000,ug))-0.85,prev);
  xv2(i)=exp(x);
  prev=x;
end


xv3=zeros(90,1);
prev=log(200);
ug=0;
mut_rate=0.05;
for i=1:length(xv2)
  r0=r0s(i)*L;
  r1=r0+(mut_rate*F);
  x=fsolve(@(x) ((1-(exp(x)-floor(exp(x))))*single_screen_power_smooth(r0,r1,floor(exp(x)),NNN,1,0.1,L,15000,ug)+ ...
                 (exp(x)-floor(exp(x)))*single_screen_power_smooth(r0,r1,floor(exp(x))+1,NNN,1,0.1,L,15000,ug))-0.85,prev);
  xv3(i)=exp(x);
  prev=x;
end

figure(1); clf;
lh(4)=loglog(r0s,xv3,'k');
hold on
lh(1)=loglog(r0s,xv,'b');
lh(2)=loglog(r0s,xv1,'g');
lh(3)=loglog(r0s,xv2,'r');
xlim([0.5e-6 20e-6]);
ylim([100 2e4]);
grid on
set(lh,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'XTick',1e-6*[0.5 1:10 20],'XTickLabel',{'    0.5','1','2','3','4','5','','','','','10', ...
                    '20'});
set(gca,'YTick',[100:100:1000 2000:1000:10000 20000],'YTickLabel',{'100','200','300','400','500','','','','',...
                    '1000','2000','3000','4000','5000','','','','','10,000','20,000'});
legend({'5%','3%','15%','10%'},'Location','NorthWest');
xlabel('Mutation rate [mut/Mb]','FontSize',18);
ylabel('Number of patients for power of 85%','FontSize',18);

% cd ~/projects/power
print_D('Number_of_samples_vs_mut_rate_5_to_1',{{'pdf'},{'png','-r180'}});



%% calc power for WGS
%%%------------------------------------------------------
clear all
close all

fnr=0.23;
F=1-fnr;
L=130;

r0=2e-6*L %used to be  6e-6*L; % changed to 10 to reply to reviewer
mut_rates=0.01:0.001:0.25;
NNN=2375714;
pv=0.1/NNN;

Ns=100;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,NNN,1,0.1,L);
  xv(i)=po;
end
xv100=xv;

Ns=300;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,NNN,1,0.1,L);
  xv(i)=po;
end
xv300=xv;

Ns=500;
xv=zeros(length(mut_rates),1);
for i=1:length(xv)
  mut_rate=mut_rates(i);
  r1=r0+(mut_rate*F);
  [po,ks]=single_screen_power_smooth(r0,r1,Ns,NNN,1,0.1,L);
  xv(i)=po;
end
xv500=xv;

figure(1); clf;
plot(mut_rates*100,xv100,'r'); hold on;
plot(mut_rates*100,xv300,'g'); hold on;
plot(mut_rates*100,xv500,'b'); hold on;

legend({'100','300','500'},'Location','SouthEast');
xlabel('% of patients with mutations','FontSize',16);
ylabel('Power','FontSize',16);
set(gca,'FontSize',14);
title(['BMR = ' num2str(1e6*r0/L) ' mut/Mb'],'FontSize',18);
grid on

p80i_97=find(xv97>=0.8,1,'first'); disp([ 97 mut_rates(p80i_97)*100]);
p80i_201=find(xv201>=0.8,1,'first'); disp([ 201 mut_rates(p80i_201)*100]);
p80i_298=find(xv298>=0.8,1,'first'); disp([ 298 mut_rates(p80i_298)*100]);


%% calc power for Lung
%%%------------------------------------------------------
clear all
close all
cd('/Users/gadgetz/Google Drive/_DROPBOX/work/matlab/gaddy/power/');

fnr=0.23;
F=1-fnr;
L=1500;

r0=12e-6*L; 
mut_rates=0.0005:0.0005:0.10;
NNN=20000;
pv=0.1/NNN;

Nss=[500 1000 2000 5000 10000 20000];
xvs={};
for ni=1:length(Nss)
    Ns=Nss(ni);
    xv=zeros(length(mut_rates),1);
    for i=1:length(xv)
        mut_rate=mut_rates(i);
        r1=r0+(mut_rate*F);
        [po,ks]=single_screen_power_smooth(r0,r1,Ns,NNN,1,0.1,L,1000);
        xv(i)=po;
    end
    xvs{ni}=xv;
end

figure(1); clf;
plot(mut_rates*100,cat(2,xvs{:})); hold on;

legend(cellstr(num2str(Nss')),'Location','SouthEast');
xlabel('% of patients with mutations','FontSize',16);
ylabel('Power','FontSize',16);
set(gca,'FontSize',14);
title(['BMR = ' num2str(1e6*r0/L) ' mut/Mb'],'FontSize',18);
grid on





