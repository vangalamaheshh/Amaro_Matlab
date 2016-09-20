function plot_segfile_horiz(segfile,lastchr)

if ~exist('lastchr','var'), lastchr = 24; end

S = load_cnseg_file(segfile);
S = reorder_struct(S,S.nprobes>=8);

len = load_chrlen;
miny = inf;
maxy = -inf;

clf

st = 0;
for c=1:lastchr, disp(c)
  Sc = reorder_struct(S,S.chr==c);
  fprintf('%d\n',slength(Sc));
  for i=1:slength(Sc)
    line(st+[Sc.start(i) Sc.end(i)],[1 1]*Sc.ratio(i),'color',[0.6 0.6 0.6],'linewidth',6);
    if Sc.ratio(i)>maxy, maxy=Sc.ratio(i); end
    if Sc.ratio(i)<miny, miny=Sc.ratio(i); end
  end
  st = st + len(c);
end

miny = miny-0.2;
maxy = maxy+0.2;

params = {'color',[0 0 0],'linewidth',1};
pos=0;
for c=1:lastchr
  pos=pos+len(c);
  line([pos pos],[miny maxy],params{:});
end
line([0 0],[miny maxy],params{:});

rectangle('position',[0 miny pos maxy-miny],'edgecolor',[0 0 0]);
ylim([miny maxy]);
xlim([0 pos]);

set(gcf,'color',[1 1 1]);
set(gca,'xtick',[],'tickdir','out');
