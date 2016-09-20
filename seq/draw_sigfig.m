function draw_sigfig(X,genes_to_show,ttypes_to_show,P);

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'margin_top',0.04);
P = impose_default_value(P,'margin_bottom',0.05);
P = impose_default_value(P,'margin_left',0.06);
P = impose_default_value(P,'margin_right',0.03);


clf,set(gca,'visible','off','position',[0 0 1 1]);

margin.t = P.margin_top; margin.b = P.margin_bottom;
margin.l = P.margin_left; margin.r = P.margin_right;
cell.h = (1-margin.t-margin.b)/length(genes_to_show);
cell.w = (1-margin.l-margin.r)/length(ttypes_to_show);
bar.h = (1/2)*cell.h;

for gi=1:length(genes_to_show), for ti=1:length(ttypes_to_show)
  xidx = find(strcmp(ttypes_to_show{ti},X.name)); gidx = find(strcmp(genes_to_show{gi},X.G{xidx}.gene));
  if isnan(xidx)||isnan(gidx), fprintf('*** WARNING: GENE+TTYPE NOT FOUND ***!\n'); continue; end
  thiscell.t = 1 - margin.t - (gi-1)*cell.h; thiscell.b = thiscell.t - cell.h;
  thiscell.l = 0 + margin.l + (ti-1)*cell.w; thiscell.r = thiscell.l + cell.w;
  if gi==1
    ttname = split(ttypes_to_show{ti},'_');
    text(thiscell.l+cell.w/2,thiscell.t+cell.h/3,ttname,'horizontalalignment','center','verticalalignment','bottom','interpreter','none');
  end
  if ti==1
    gname = genes_to_show{gi};
    text(thiscell.l-cell.w/10,thiscell.b+cell.h/2,gname,'horizontalalignment','right','interpreter','none');
  end
  s = max(0,min(6,-log10([X.G{xidx}.pmaxCV(gidx) X.G{xidx}.pcons(gidx) X.G{xidx}.pclust(gidx) X.G{xidx}.pmax(gidx)])));
  eff = [X.G{xidx}.effCV(gidx) X.G{xidx}.effcons(gidx) X.G{xidx}.effclust(gidx)];
  effmax = [4 2.5 10]; if ~grepm('pancan|leave|best|union',ttypes_to_show(ti)), effmax(1) = 20; end
  fraceff = max(0,min(1,eff./effmax));
  % cell background color: greyscale proportional to %patients, capped at white=0%, black=10%
  pctpat = 100*X.G{xidx}.npat(gidx)/X.npat(xidx); pctpat=min(10,max(0,pctpat)); bkgdcolor=1-pctpat/10*ones(1,3);
  rectangle('position',[thiscell.l thiscell.b cell.w cell.h],'facecolor',bkgdcolor,'edgecolor',0.9*bkgdcolor+0.1*[0.5 0.5 0.5]);
  % bar: total length = -log10(pmax), capped at 6
  thisbar.width = cell.w*(s(4)/6); thisbar.b = thiscell.b+(bar.h/2);
  if thisbar.width>0, rectangle('position',[thiscell.l thisbar.b thisbar.width bar.h],'facecolor',[0.5 0.5 0.5]); end
  % bar segments: fractional lengths = -log10(pcv,pcons,pclust), each capped at 6
  % bar colors: R, G, B, with saturation ~= effect size
  fracwidth = s(1:3)/sum(s(1:3)); thisbar.segwidth = thisbar.width * fracwidth;
  x1 = thiscell.l; x2 = x1 + thisbar.segwidth(1); x3 = x2 + thisbar.segwidth(2);
  if thisbar.segwidth(1)>0, rectangle('position',[x1 thisbar.b thisbar.segwidth(1) bar.h],'facecolor',[1 1-fraceff(1) 1-fraceff(1)]); end
  if thisbar.segwidth(2)>0, rectangle('position',[x2 thisbar.b thisbar.segwidth(2) bar.h],'facecolor',[1-fraceff(2) 1 1-fraceff(2)]); end
  if thisbar.segwidth(3)>0, rectangle('position',[x3 thisbar.b thisbar.segwidth(3) bar.h],'facecolor',[1-fraceff(3) 1-fraceff(3) 1]); end
end,end




