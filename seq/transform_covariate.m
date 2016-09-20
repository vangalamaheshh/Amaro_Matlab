function x2 = transform_covariate(G,xfld,yfld)
% xfld  e.g. 'hiC'   == covariate
% yfld  e.g 'rate'   == output we're trying to fit/model

smoothfrac = 0.1 ;  % results depend on this.  range of 0.2-0.5 seems to give decent results.  lower is noisier.

%G = reorder_struct(G,G.Nbkgd>1e7);

demand_fields(G,{'name';xfld;yfld});

[q idx1] = reorder_struct(G,~isnan(getfield(G,xfld))); [q idx2] = sort_struct(q,xfld);
rate = getfield(q,yfld); lograte = log10(rate); n = length(rate); cv = getfield(q,xfld);
x = cv;
y = rate;

clf,subplot('position',[0.05 0.2 0.4 0.7]),hold on
scatter(x,y,10,[0.7 0.7 0.7]);
g = grepm('^OR\d',q.name);
scatter(x(g),y(g),10,[1 0 0],'filled');
xl=xlim;yl=ylim; for yi=0:0.5:7, line(xl,yi*[1 1],'color',[0.3 0.3 0.3]); end; ylim([0 7]);
xadj = diff(xl)*0.01;
% smooth
smoothrate = smooth(rate,smoothfrac*n,'loess');
% display smoothed curve
plot(x,smoothrate,'-','color',[0 0 0]);
% highlight some genes
genes = {'CSMD3','CSMD1','RYR2','RYR3','TPTE','PCLO','TTN'};
g = ismember(q.name,genes);
scatter(x(g),y(g),30,[0 0 1],'filled');
textfit(x(g)+xadj,y(g),q.name(g),'color',[0 0 1],'fontsize',8);
genes = {'SF3B1','MYD88','HRAS','KRAS','NRAS','TP53'};
g = ismember(q.name,genes);
scatter(x(g),y(g),30,[0 0.7 0],'filled');
textfit(x(g)+xadj,y(g),q.name(g),'color',[0 0.7 0],'fontsize',8);

hold off,fp,title(xfld,'interpreter','none');

% transform
newfld = [xfld '_transformed']; x2 = nan(slength(G),1); x2(idx1(idx2)) = smoothrate;G = setfield(G,newfld,x2);

% measure of how well it's working
discrim_ratio = max(x2)/min(x2);
fprintf('discrim_ratio = %f\n',discrim_ratio);


[q idx1] = reorder_struct(G,~isnan(getfield(G,newfld))); [q idx2] = sort_struct(q,newfld);
rate = getfield(q,yfld); lograte = log10(rate); n = length(rate); cv = getfield(q,newfld);
x = cv;
y = rate;

subplot('position',[0.55 0.2 0.4 0.7]),hold on
scatter(x,y,10,[0.7 0.7 0.7]);
g = grepm('^OR\d',q.name);
scatter(x(g),y(g),10,[1 0 0],'filled');
xl=xlim;yl=ylim; for yi=0:0.5:7, line(xl,yi*[1 1],'color',[0.3 0.3 0.3]); end; ylim([0 7]);
xadj = diff(xl)*0.01;
% smooth
smoothrate = smooth(rate,smoothfrac*n,'loess');
% display smoothed curve
plot(x,smoothrate,'-','color',[0 0 0]);
% highlight some genes
genes = {'CSMD3','CSMD1','RYR2','RYR3','TPTE','PCLO','TTN'};
g = ismember(q.name,genes);
scatter(x(g),y(g),30,[0 0 1],'filled');
textfit(x(g)+xadj,y(g),q.name(g),'color',[0 0 1],'fontsize',8);
genes = {'SF3B1','MYD88','HRAS','KRAS','NRAS','TP53'};
g = ismember(q.name,genes);
scatter(x(g),y(g),30,[0 0.7 0],'filled');
textfit(x(g)+xadj,y(g),q.name(g),'color',[0 0.7 0],'fontsize',8);

hold off,fp,title(newfld,'interpreter','none');


set(gcf,'position',[384   312   751   380])






