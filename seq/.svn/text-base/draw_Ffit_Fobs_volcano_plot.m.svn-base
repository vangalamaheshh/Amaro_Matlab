function draw_Ffit_Fobs_volcano_plot(G,P)

if ~exist('P','var'), P=[];end
P = impose_default_value(P,'volcano_left_boundary',-0.5);
P = impose_default_value(P,'volcano_top_boundary',18);
P = impose_default_value(P,'mark_bad_genes', ...
   {'TTN','MUC16','LRP1B','CSMD1','CSMD3','RYR2','RYR3','PCLO','TPTE','OR*'});
P = impose_default_value(P,'mark_good_genes', ...
   {'TP53','KRAS','NRAS','SF3B1','MYD88','EGFR','ERBB2','NF1','RB1','APC','ATM','PIK3CA','VHL','PTEN','STK11','NFE2L2','KEAP1','CDKN2A'});

demand_fields(G,{'name','Fobs','Ffit','p'});
G = make_numeric(G,{'Fobs','Ffit','p'});
ng = slength(G);

G.Fratio = G.Fobs./G.Ffit;
G.Fratio(G.Ffit==0) = 100;
G.logFratio = log10(G.Fratio);
G.logFratio(G.Fratio<=0) = P.volcano_left_boundary;
G.logFratio(G.logFratio<P.volcano_left_boundary) = P.volcano_left_boundary;

G.score = -log10(G.p); G.score(G.score>P.volcano_top_boundary) = P.volcano_top_boundary;

x = G.logFratio;
y = G.score;

gene_marker_sizes = repmat(20,ng,1);
gene_marker_colors = repmat([0.7 0.7 0.7],ng,1);

figure(1);clf;
scatter(x,y,gene_marker_sizes,gene_marker_colors);
hold on
xlabel('log10 Fobs/Ffit','fontsize',20);
ylabel('score = -log10 p\_MutSig','fontsize',20);

idx = find((y>=2 & x>=1) | (y>=4 & x>=0.8) | (y>=10 & x>=0.6));
xl=xlim; dx=0.08*(diff(xl)/10);
textfit(x(idx)+dx,y(idx),G.name(idx))
return

%%%%% mark some genes

genes_to_mark = [as_column(P.mark_bad_genes);as_column(P.mark_good_genes)];

if isfield(G,'q')  % if MutSig q-values are available, mark all significant genes
  G = make_numeric(G,'q');
  sig_genes = G.name(G.q<=P.mutsig_q_cutoff);
  gidx = listmap(sig_genes,G.name);
  scatter(x(gidx),y(gidx),20,[0.5 0.5 0.5],'filled');
end

if isfield(G,'rank')   % if MutSig rankings are available, mark the top N genes
  G = make_numeric(G,'rank');
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
if 0 && isfield(G,'q')  % if MutSig q-values are available, mark all significant genes
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
scatter(x(gidx),y(gidx),gene_marker_sizes,gene_marker_colors,'filled');
xl=xlim; dx=0.08*(diff(xl)/10);
textfit(x(gidx)+dx,y(gidx),G.name(gidx),'color',gene_font_colors,'fontsize',gene_font_sizes);

hold off,set(gcf,'color',[1 1 1],'position',[296 37 904 786]);

