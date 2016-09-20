function create_gbrowse_qqplots(P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'runfile','*required*');
P = impose_default_value(P,'datasetdir','*required*');

fprintf('Loading runfile.\n');
load(P.runfile,'X');
for i=1:slength(X)
  if ~isfield(X.G{i},'pmid') && isfield(X.G{i},'p'), X.G{i} = rename_field(X.G{i},'p','pmid'); end
  if ~isfield(X.G{i},'pmidCV') && isfield(X.G{i},'pCV'), X.G{i} = rename_field(X.G{i},'pCV','pmidCV'); end
end

fprintf('Creating QQ plots... ');
outdir = [P.datasetdir '/qqplots']; ede(outdir);
figure(1),for i=1:slength(X),fprintf('%d/%d ',i,slength(X));
  if strcmp(X.name{i},'union') || strcmp(X.name{i},'best_ttype')
    exclude = 'union|best_ttype|leaveout|\_|COAD|READ'; % don't count leaveouts or subtypes
    if strcmp(X.name{i},'best_ttype'), exclude = ['pan|' exclude]; end  % or pancan (if this isn't union)
    xi = grepvi(exclude,X.name,1);
    pmid = cell(length(xi),1); for j=1:length(xi), pmid{j} = X.G{xi(j)}.pmid; end
    pmid = cat(1,pmid{:});
    vals_for_plot = {pmid};
    txts_for_plot = {'joint3'};
    cols_for_plot = {[1 0 1]};
    tots_for_plot = [sum(X.G{i}.q<=0.1)];
  else % one of the original runs
    vals_for_plot = {X.G{i}.pclust,X.G{i}.pcons,X.G{i}.pjoint,X.G{i}.pmidCV,X.G{i}.pmid};
    txts_for_plot = {'clust','cons','joint2','cv','joint3'};
    cols_for_plot = {[0 0 1],[0 0.7 0],[1 0 0],[0 1 1],[1 0 1]};
    pmax = {X.G{i}.pclust,X.G{i}.pcons,X.G{i}.pjoint,X.G{i}.pmaxCV,X.G{i}.pmax};
    tots_for_plot = nan(length(pmax),1);
    for j=1:length(pmax)
      q = calc_fdr_value(pmax{j});
      tots_for_plot(j) = sum(q<=0.1);
    end
  end
  % fill in missing values with random numbers
  for j=1:length(vals_for_plot), idx = find(isnan(vals_for_plot{j})); vals_for_plot{j}(idx) = rand(length(idx),1); end
  clf,qq(vals_for_plot{:}); xlim([0 5]); ylim([0 7]); yl=ylim;
  text(0.1,yl(2)-0.75,X.name{i},'fontsize',60,'color',[0 0 0],'interpreter','none');
  for j=1:length(vals_for_plot)
    txt = sprintf('%s (%d)',txts_for_plot{j},tots_for_plot(j));
    text(0.1,yl(2)-0.85-0.50*j,txt,'fontsize',30,'color',cols_for_plot{j});
  end
  set(gca,'position',[0.07 0.05 0.90 0.92]); set(gcf,'color',[1 1 1]);
  print('-dpng','-r40',[outdir '/ttype' num2str(i) '.qq.png']); % lowres for webpage
end,fprintf('\n');






