function M = compute_nonsilent_mle_bounded_by_silent_ci(M,g,ci)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'silent_ci_to_use',ci);

%fprintf('compute_nonsilent_mle_bounded_by_silent_ci... ');

%%% for each gene:
% (1) get posterior distribution on silent rate
% (2) find 95% confidence interval
% (3) find MLE for nonsilent rate
% (4) if (3) is within (2), then use (3)
% (4) else, bound (2) within (3)


%%  first approximation: ignore categories

mutcats = 1:M.TOT-M.NUM_INDEL_CLASSES;
covcats = M.TOT;

G = [];
G.name = M.gene.name;
G.nsil = sum(sum(M.n_silent(:,mutcats,:),3),2);
G.Nsil = round(sum(sum(M.N_sil_cov(:,covcats,:),3),2));

G.nnon = sum(sum(M.n_nonsilent_ignoring_null_categ(:,mutcats,:),3),2);
G.Nnon = round(sum(sum(M.N_non_cov(:,covcats,:),3),2));

alpha = 1-P.silent_ci_to_use;
G.mle_non = 1e6*binofit(G.nnon,G.Nnon,alpha);
[tmp1 tmp2] = binofit(G.nsil,G.Nsil,alpha);
G.ci_sil = 1e6*tmp2;



G.mu_non = G.mle_non;
idx = find(G.ci_sil(:,2) < G.mle_non);
G.mu_non(idx) = G.ci_sil(idx,2);

pval = 1-binocdf(G.nnon(g)-1,G.Nnon(g),G.mu_non(g)/1e6);


%%  graphical exploration

mu_test_range = 0:0.5:2000;
post_sil = binopdf(G.nsil(g),G.Nsil(g),mu_test_range/1e6);
post_non = binopdf(G.nnon(g),G.Nnon(g),mu_test_range/1e6);
post_sil_norm = post_sil/max(post_sil);
post_non_norm = post_non/max(post_non);
crop = find(post_sil_norm>0.01 | post_non_norm>0.01,1,'last');
post_sil_norm = post_sil_norm(1:crop); post_non_norm = post_non_norm(1:crop); mu_range = mu_test_range(1:crop);

clf,hold on
plot(mu_range,post_sil_norm,'-','color',[0 0.6 0]);
plot(mu_range,post_non_norm,'-','color',[1 0 0]);
yt = 1.5; ylim([0 yt]); ylabel('post/max(post)');
xlim(mu_range([1 end])); xlabel('mutation rate (/Mb)');
set(gca,'tickdir','out');
title(G.name{g},'fontsize',20);
rgt = {'horizontalalignment','right'}; ctr = {'horizontalalignment','center'};
tx = mu_range(end)*0.95; ty = yt-0.2;
x1=G.mle_non(g);line(x1*[1 1],[0 ty-0.18],'color',[1 0 0]); text(x1,ty-0.16,'MLE',ctr{:},'color',[1 0 0]);
x2=G.ci_sil(g,2);line(x2*[1 1],[0 ty-0.25],'color',[0 0.6 0]);text(x2,ty-0.23,sprintf('%2.0f%%ci(hi)',100*(1-alpha)),ctr{:},'color',[0 0.6 0]);
x=min(x1,x2);line(x*[1 1],ty-[0.11 0.03],'color',[0 0 0]); text(x,ty,'mu(g)',ctr{:},'color',[0 0 0]);
line(x*[1 1.02],ty-[0.11 0.08],'color',[0 0 0]);
line(x*[0.98 1],ty-[0.08 0.11],'color',[0 0 0]);
text(tx,yt-0.05,sprintf('silent: n = %d  N = %d',G.nsil(g),G.Nsil(g)),rgt{:});
text(tx,yt-0.11,sprintf('nonsilent: n = %d  N = %d',G.nnon(g),G.Nnon(g)),rgt{:});
if pval<0.001, fmt = '%d'; else fmt = '%f'; end
text(tx,yt-0.17,sprintf(['p = ' fmt],pval),rgt{:});
hold off



if nargout==0, clear M, end
return





% compute adjusted BMRs



M.mutrate.per_gene_BMR_correction = ones(M.ng,1);

global_non_bmr = M.mutrate.tot.hat;

%idx = find(
