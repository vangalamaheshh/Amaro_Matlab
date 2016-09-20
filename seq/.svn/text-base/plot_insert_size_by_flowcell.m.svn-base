function plot_insert_size_by_flowcell(T,N,flowcell_order,randseed,yl,individual_name)
% Mike Lawrence 2009-06-25
%      modified 2009-10-15

if ~exist('randseed','var') | isempty(randseed), randseed = 1234; end
if ~exist('yl','var'), yl=[]; end

if ~exist('individual_name','var')
  ttext = [];
else
  ttext = [individual_name ' '];
end

T.fc = regexprep(T.lanes.PU,'AAXX.*','');
N.fc = regexprep(N.lanes.PU,'AAXX.*','');

if exist('flowcell_order','var') & ~isempty(flowcell_order)
  ufc = flowcell_order;
  needed = unique([T.fc(T.good);N.fc(N.good)]);
  missing = setdiff(needed,ufc);
  if ~isempty(missing)
    fprintf('"flowcell_order" is missing the following entries:');
    for i=1:length(missing), fprintf(' %s',missing{i}); end
    fprintf('\n');
    error('Please add them!');
  end
else
  ufc = unique([T.fc;N.fc]);
end

T.fc_no = listmap(T.fc,ufc);
N.fc_no = listmap(N.fc,ufc);

legtext = {'tumor lanes';'normal lanes'};
legcol = [0 0 0; 0.6 0.6 0.6];
tcol = zeros(length(T.good),3);
ncol = 0.6*ones(length(N.good),3);

if isfield(T,'outliers') && ~isempty(T.outliers)
  legtext = [legtext;'outlier tumor lanes'];  legcol = [legcol; 1 0 0];
  tcol(T.outliers,1) = 1; % red
end
if isfield(N,'outliers') && ~isempty(N.outliers)
  legtext = [legtext;'outlier normal lanes'];  legcol = [legcol; 1 0.6 0.6];
  ncol(N.outliers,1) = 1; % pink
end
if isfield(T,'blacklisted') && ~isempty(T.blacklisted)
  legtext = [legtext;'blacklisted tumor lanes'];  legcol = [legcol; 0 0 1];
  tcol(T.blacklisted,[1 2]) = 0;
  tcol(T.blacklisted,3) = 1; % blue
end
if isfield(N,'blacklisted') && ~isempty(N.blacklisted)
  legtext = [legtext;'blacklisted normal lanes'];  legcol = [legcol; 0.5 0.5 1];
  ncol(N.blacklisted,[1 2]) = 0.5;
  ncol(N.blacklisted,3) = 1; % light blue
end




figure(1);clf
catplot([T.adjmean(T.good);N.adjmean(N.good)] ,...
        [T.fc_no(T.good);N.fc_no(N.good)] ,...
        [tcol;ncol],ufc,0.5,'o',randseed,yl);

title([ttext 'insert size'],'fontsize',20,'interpreter','none');
ylabel('adjusted mean','fontsize',20);

set(gca,'position',[0.13 0.16 0.78 0.73]);
xl = xlim; yl = ylim;
text(xl(1)+xl(2)/2-0.1,yl(1)-0.12*(yl(2)-yl(1)),'flowcell','fontsize',20,...
   'horizontalalign','center','verticalalign','top');

if isempty(T.good) && isempty(N.good)
  text(0.5,0.5,'No paired reads','fontsize',20,'horizontalalignment','center');
else
  l=legend(legtext{:},'Location','best');
  a = get(l,'children');
  idx = size(legcol,1);
  for i=1:3:length(a)
    set(a(i),'markersize',6,'markeredgecolor',legcol(idx,:),'markerfacecolor',legcol(idx,:));
    idx=idx-1;
  end
end


