function fa=make_pharma_figure(D,C,gsupids,regs,pvs,col)

fa=NaN;
vis='off';

figure;
fa=gcf;

if ~isempty(regs)
  k=1;
  sp=plot_snp_score('',C,pvs,[],-0.25,0,0);
  sp=sp{k};
  regs=regs{k};
  pvs=pvs{k};
else
  figure;
  axis([0 length(pvs)+1 0 max(pvs)]);
  w=0.1;
  for j=1:length(pvs)
    bh(j)=patch([j-w/2 j+w/2 j+w/2 j-w/2 j-w/2],[0 0 pvs(j) pvs(j) 0],[0 0 1]);
  end
%  bar(1:length(pvs),pvs,0.1);
  sp=gca;
end

axsp=[get(sp,'xlim')' get(sp,'ylim')' get(sp,'zlim')'];
axes(sp);
ax=axis;
pan=patch(ax([ 1 2 2 1 1]), ax([3 3 4 4 3]), [ 0 0 0 0 0]);
set(pan,'FaceColor',[0 0 0],'FaceAlpha',0.1);
ftmp=gcf;


figure(fa);
clf;
%sz=[length(pvs) 20 length(pvs)];
sz=[100 20 100];
axis([ 0 sz(1) 0 sz(2) 0 sz(3)]);
ax=axis;

pt={};
pt{1}=[ 0 0 0; sz(1) 0 0; sz(1) sz(2) 0; 0 sz(2) 0; 0 0 0];

for i=1:length(pt)
  ph(i)=patch(pt{i}(:,1),pt{i}(:,2),pt{i}(:,3),i);
end
set(ph,'FaceAlpha',0.1);
set(ph,'Visible',vis);

set(gca,'CameraPosition',[sz(1)*10 sz(2)*50 -sz(3)*10]);
set(gca,'CameraTarget',[mean(ax(1:2)) mean(ax(3:4)) mean(ax(5:6))]);
set(gca,'CameraUpVector',[0 1 0]);

put_in_3d(get(sp,'Children'),pt{1}); 
close(ftmp);

for i=1:length(gsupids)
  disp(i);
  s=D.gsupdat(gsupids(i),:);
  
  figure(ftmp); clf;
  mins=min(s);
  maxs=max(s);
  delta=0.1;
  x=mins:delta:maxs;
  hst=histc(s,x);
  %    bh=bar(x,hst);
  bh=[];
  for j=1:length(hst)
    bh(j)=patch([ x(j)-delta/2 x(j)+delta/2 x(j)+delta/2 x(j)-delta/2 x(j)-delta/2],[0 0 log10(hst(j)+1) log10(hst(j)+1) 0],col);
  end
  set(bh,'FaceAlpha',0.1,'LineStyle','none');
  min(-4,mins-delta/2)
  t1=min(-4,mins-delta/2);
  t2=max(4,maxs+delta/2);
  
%  lh1=line([t1 4],[0 0],'Color',[0.8 0.8 0.8]);
  lh=patch([t1 t2 t2 t1 t1],[0 0 -0.05 -0.05 0],col); set(lh,'FaceAlpha',0.1,'LineStyle','none');
  
  sigj=find(x<-4);
  if ~isempty(sigj)
    set(bh(sigj),'FaceColor',col,'FaceAlpha',0.8);
  end
  ax=axis;
  axis([-4 4 0 ax(4)]);
%  ax=axis;
%  lh=patch(ax([ 1 2 2 1 1]),ax([ 3 3 4 4 3]),[0 1 0]); set(lh,'FaceAlpha',0.2,'LineStyle','none');
  ta=gca;
  
  figure(fa);
  if ~isempty(regs)
    scx=sz(1)/length(pvs);
    pt{i+1}=[ regs(i).peak*scx 0 sz(3); regs(i).peak*scx 0 0; regs(i).peak*scx sz(2) 0; ...
              regs(i).peak*scx sz(2) sz(3); regs(i).peak*scx 0 sz(3); ];
  else
    scx=sz(1)/(length(pvs)+1);    
    pt{i+1}=[ i*scx 0 sz(3); i*scx 0 0; i*scx sz(2) 0; ...
              i*scx sz(2) sz(3); i*scx 0 sz(3); ];
  end
  
  ph(i+1)=patch(pt{i}(:,1),pt{i}(:,2),pt{i}(:,3),i+1);
  set(ph(i+1),'FaceAlpha',0.1);
  set(ph(i+1),'Visible',vis);
  
  put_in_3d(get(ta,'Children'),pt{i+1});
%  keyboard
  
  close(ftmp);
  %    pause
end

