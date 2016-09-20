function plot_TN_isz_distribs(T,N,individual_name,P)
% Mike Lawrence 2009

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'xmin',0);
P = impose_default_value(P,'xmax',500);

xrange = [P.xmin P.xmax];

if ~exist('individual_name','var')
  ttext = [];
else
  ttext = [individual_name ' '];
end

% TUMOR LANES

subplot(2,1,1);hold on

legtext = {'tumor lanes'}; legcol = [0 0 0];
plot(T.datnorm(:,T.good),'color',[0 0 0]);

if isfield(T,'outliers') && ~isempty(T.outliers)
  plot(T.datnorm(:,T.outliers),'color',[1 0 0]);
  legtext = [legtext; 'outlier tumor lanes']; legcol = [legcol; 1 0 0];
end
if isfield(T,'blacklisted') && ~isempty(T.blacklisted)
  plot(T.datnorm(:,T.blacklisted),'color',[0 0 1]);
  legtext = [legtext; 'blacklisted tumor lanes']; legcol = [legcol; 0 0 1];
end
xlim(xrange);xlabel('insert size','fontsize',15);title([ttext 'tumor lanes'],'fontsize',20,'interpreter','none');
set(gca,'ytick',[]);ylabel('fraction of lane','fontsize',15);
if isempty(T.good), text((P.xmin+P.xmax)/2,0.5,'No paired reads','fontsize',20,'horizontalalignment','center'); end

l=legend(legtext{:},'Location','northwest');
a = get(l,'children');
idx = size(legcol,1);
for i=2:3:length(a)
  set(a(i),'Color',legcol(idx,:));
  idx=idx-1;
end

hold off;


% NORMAL LANES

subplot(2,1,2);hold on

legtext = {'normal lanes'}; legcol = [0.6 0.6 0.6];
plot(N.datnorm(:,N.good),'color',[0.6 0.6 0.6]);

if isfield(N,'outliers') && ~isempty(N.outliers)
  plot(N.datnorm(:,N.outliers),'color',[1 0.6 0.6]);
  legtext = [legtext; 'outlier normal lanes']; legcol = [legcol; 1 0.6 0.6];
end
if isfield(N,'blacklisted') && ~isempty(N.blacklisted)
  plot(N.datnorm(:,N.blacklisted),'color',[0.5 0.5 1]);
  legtext = [legtext; 'blacklisted normal lanes']; legcol = [legcol; 0.5 0.5 1];
end
xlim(xrange);xlabel('insert size','fontsize',15);title([ttext 'normal lanes'],'fontsize',20,'interpreter','none');
set(gca,'ytick',[]);ylabel('fraction of lane','fontsize',15);
if isempty(N.good), text((P.xmin+P.xmax)/2,0.5,'No paired reads','fontsize',20,'horizontalalignment','center'); end

l=legend(legtext{:},'Location','northwest');
a = get(l,'children');
idx = size(legcol,1);
for i=2:3:length(a)
  set(a(i),'Color',legcol(idx,:));
  idx=idx-1;
end

hold off
