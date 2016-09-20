function scatter_width_vs_mean(T,N,individual_name)
% Mike Lawrence 2009-06-25
%      modified 2009-10-15

if ~exist('individual_name','var')
  ttext = [];
else
  ttext = [individual_name ' '];
end

figure(1),clf,hold on
legtext = {};
legtext = [legtext;'tumor lanes'];
scatter(T.adjmean(T.good),T.width(T.good),40,[0 0 0],'o','filled');
legtext = [legtext;'normal lanes'];
scatter(N.adjmean(N.good),N.width(N.good),50,[0.6 0.6 0.6],'^','filled');
if isfield(T,'outliers') && ~isempty(T.outliers)
  scatter(T.adjmean(T.outliers),T.width(T.outliers),40,[1 0 0],'o','filled');
  legtext = [legtext;'outlier tumor lanes'];
end
if isfield(N,'outliers') && ~isempty(N.outliers)
  scatter(N.adjmean(N.outliers),N.width(N.outliers),50,[1 0.6 0.6],'^','filled');
  legtext = [legtext;'outlier normal lanes'];
end
if isfield(T,'blacklisted') && ~isempty(T.blacklisted)
  scatter(T.adjmean(T.blacklisted),T.width(T.blacklisted),40,[0 0 1],'o','filled');
  legtext = [legtext;'blacklisted tumor lanes'];
end
if isfield(N,'blacklisted') && ~isempty(N.blacklisted)
  scatter(N.adjmean(N.blacklisted),N.width(N.blacklisted),50,[0.5 0.5 1],'^','filled');
  legtext = [legtext;'blacklisted normal lanes'];
end



title([ttext 'insert size'],'fontsize',20,'interpreter','none');
xlabel('adjusted mean','fontsize',20); ylabel('peak width','fontsize',20);

if isempty(T.good) && isempty(N.good)
  text(0.5,0.5,'No paired reads','fontsize',20,'horizontalalignment','center');
else
  l=legend(legtext{:},'Location','best');
  a = get(l,'children');
  for i=1:length(a)
    b = get(a(i),'children');
    set(b,'markersize',6);
  end
end

hold off
