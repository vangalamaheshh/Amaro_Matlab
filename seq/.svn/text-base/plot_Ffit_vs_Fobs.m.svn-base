function plot_Ffit_vs_Fobs(G,P)

if ~exist('P','var'), P=[];end
P = impose_default_value(P,'axismax',10);
P = impose_default_value(P,'report_top_n_mutsig_genes',40);
P = impose_default_value(P,'mutsig_q_cutoff',0.1);
P = impose_default_value(P,'min_Fdiff',4);
P = impose_default_value(P,'mark_bad_genes', ...
   {'TTN','MUC16','LRP1B','CSMD1','CSMD3','RYR2','RYR3','PCLO','TPTE'});
P = impose_default_value(P,'mark_good_genes', ...
   {'TP53','KRAS','NRAS','SF3B1','MYD88','EGFR','ERBB2','NF1','RB1','APC','ATM','PIK3CA','VHL','PTEN','STK11','NFE2L2','KEAP1','CDKN2A'});

demand_fields(G,{'name','Fobs','Ffit'});
G = make_numeric(G,{'Fobs','Ffit'});
name=G.name;
x=G.Fobs;
y=G.Ffit;

idx=find(x>P.axismax);
for j=1:length(idx),i=idx(j);name{i}=sprintf('%s (%.0f)',name{i},x(i));end;
x(x>P.axismax)=P.axismax;y(y>P.axismax)=P.axismax;

figure(1);clf
scatter(x,y,20,[0.7 0.7 0.7]);
hold on
xl=xlim;yl=ylim;mx=min(max([xl(2),yl(2)]),P.axismax);xlim([0 mx]);ylim([0 mx]);line([0 mx],[0 mx],'color',[0 0 0]);
xlabel('F observed','fontsize',20); ylabel('F predicted','fontsize',20);

%%%%% mark some genes

if isfield(G,'q') && isfield(G,'rank')  % if MutSig q-values are available, mark significant genes
  G = make_numeric(G,'q');
  sig_genes = G.name(G.q<=P.mutsig_q_cutoff);
  gidx = listmap(sig_genes,G.name);
  scatter(x(gidx),y(gidx),20,[0.5 0.5 0.5],'filled');
  G = make_numeric(G,'rank');  % also label the top N genes
  genes_to_mark = G.name(G.rank<=P.report_top_n_mutsig_genes);
else   % otherwise just show the complete canned bad/good lists
  genes_to_mark = [as_column(P.mark_bad_genes);as_column(P.mark_good_genes)];
end

gene_marker_colors = repmat([0.5 0.5 0.5],length(genes_to_mark),1);
gene_marker_sizes = 50*ones(length(genes_to_mark),1);
gene_font_colors = repmat([0 0 0],length(genes_to_mark),1);
gene_font_sizes = 10*ones(length(genes_to_mark),1);

% mark "bad" genes red
idx = find(ismember(genes_to_mark,P.mark_bad_genes));
gene_marker_colors(idx,:) = repmat([1 0 0],length(idx),1);
gene_font_colors(idx,:) = repmat([0.6 0 0],length(idx),1);
gene_font_sizes(idx) = 12;

% mark "good" genes green
idx = find(ismember(genes_to_mark,P.mark_good_genes));
gene_marker_colors(idx,:) = repmat([0 0.7 0],length(idx),1);
gene_font_colors(idx,:) = repmat([0 0.3 0],length(idx),1);
gene_font_sizes(idx) = 12;

% write names of significant genes of minimum Fdiff=(Fobs-Ffit) that didn't make the top N
if isfield(G,'q')  % if MutSig q-values are available, mark all significant genes
  G = make_numeric(G,'q');
  sig_genes = G.name(G.q<=P.mutsig_q_cutoff);
  G.Fdiff = G.Fobs-G.Ffit;
  f_genes = G.name(G.Fdiff>=P.min_Fdiff);
  sigf_genes = intersect(sig_genes,f_genes);
  missed_sigf_genes = setdiff(sigf_genes,genes_to_mark);
  genes_to_mark = [genes_to_mark; missed_sigf_genes];
  gene_marker_colors = [gene_marker_colors; repmat([0.5 0.5 0.5],length(missed_sigf_genes),1)];
  gene_marker_sizes = [gene_marker_sizes; repmat(20,length(missed_sigf_genes),1)];
  gene_font_colors = [gene_font_colors; repmat([0 0 0],length(missed_sigf_genes),1)];
  gene_font_sizes = [gene_font_sizes; repmat(10,length(missed_sigf_genes),1)];
end

gidx = listmap(genes_to_mark,G.name);
gidx(isnan(gidx))=[];
scatter(x(gidx),y(gidx),gene_marker_sizes,gene_marker_colors,'filled');
textfit(x(gidx)+0.08,y(gidx),name(gidx),'color',gene_font_colors,'fontsize',gene_font_sizes);

hold off,set(gcf,'color',[1 1 1],'position',[296 37 904 786]);

