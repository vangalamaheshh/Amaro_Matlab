function ph = catplot(v,catg,cols,labels,addrnd,shape,randseed,yl,P)
% catplot(v,catg,cols,labels,addrnd,shape,randseed,yl)
%
% v = y-axis position
% catg = category number
% cols = colors
% labels = name of each category
% addrnd = how much x-randomness to add (e.g. 0.5)
% shape = plot shape e.g. 'o'
% randseed = random seed (if nan, points are plotted sequentially)
% yl = ylim
%
% returns handles to all points drawn
% 
% Gaddy Getz 2008

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'fontsize', 10);

if exist('randseed','var')
  if ~isnan(randseed)
    rand('twister',randseed);
    verbose(['setting randseed to ' num2str(randseed)],30);
  end
else
  rand('twister',5489);
end

wd=0.02;
hgt=0.04;
alpha=0.5;
if ~exist('addrnd','var')
  addrnd=0;
end

[u,ui,uj]=unique(catg);
ph=cell(1,length(u));
plot(v);
enlarge_axis([0 0.5]);
ax=axis;
axyt=get(gca,'YTick');
axytl=get(gca,'YTickLabel');
axpos = get(gca,'position');
clf;
set(gca,'position',axpos);
u(isnan(u))=[];
for i=1:length(u)
  x=i;
  ini=find(uj==i);
  p=[];
  for j=1:length(ini)
    k = ini(j);
    y=v(k);

    %%%%   FEATURE ADDED  2009-07-14 Mike Lawrence
    %%%%   if randseed is "nan", then the dots are spaced evenly.

    if ~isnan(randseed)
      r = (rand(1,1)-0.5)*addrnd;
    else
      r = (-0.5 + (j/(length(ini)+1))) * addrnd;
    end

%%%%   BUGFIX  2009-06-24 Mike Lawrence
%%%%   --> was using "i" as index to cols;
%%%%       should be using v(k) instead!

    if strmatch(shape,'patch')
      p(j)=patch([ x-wd+r x+wd+r x+wd+r x-wd+r x-wd+r],[y-hgt y-hgt y+hgt y+hgt y-hgt],cols(k,:));
      set(p(j),'FaceAlpha',alpha,'EdgeColor','none');
    else
      p(j)=plot(x+r,y,shape);
      set(p(j),'MarkerFaceColor',cols(k,:),'MarkerEdgeColor',cols(k,:));
      hold on;
    end
  end
  ph{k}=p;
end

if ~isempty(v) && ~isempty(catg)
  axis([ 0.5 length(u)+0.5 ax(3:4)] );
  set(gca,'YTick',axyt,'YTickLabel',axytl,'fontsize',P.fontsize);
end

if exist('yl','var') && ~isempty(yl), ylim(yl); end

%%%%   BUGFIX  2009-06-24 Mike Lawrence
%%%%   --> was using labels, unordered.
%%%%       should be using labels(catg(ui)) instead.

tick = 1:length(u);
label = labels(catg(ui));
if length(u)<10
  set(gca,'XTick',tick,'XTickLabel',label,'fontsize',P.fontsize);
else
  xticklabel_rotate_simplified(tick,90,label,'verticalalign','middle','fontsize',P.fontsize);
end
  



