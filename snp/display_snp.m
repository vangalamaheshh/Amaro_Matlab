function h=display_snp(CL21,regs)

figure(1); clf;
if isfield(CL21,'raw');
  h=subplot(1,3,1);
  x=CL21.raw;
  if abs(mean(x(:)))>0.1
    x(x<0.1)=0.1;
    x=log2(x)-1;
  end
  imagesc(x);
  ax(1)=gca;
  bluepink
  hc(1)=colorbar;
  zoom yon
else
  h=subplot(1,2,1);
  imagesc(CL21.dat);
  ax(1)=gca;
  bluepink
  hc(1)=colorbar;
  zoom yon
  caxis([-1 1]);
end

if isfield(CL21,'sm3')
  subplot(1,3,2);
  imagesc(CL21.sm3);
  ax(2)=gca;
  bluepink
  hc(2)=colorbar;
  zoom yon
  caxis([-1 1]);
else
  ax(2)=NaN;
end

for i=1:max(CL21.chrn)
  chrnpos(i)=round(mean(find(CL21.chrn==i)));
end

if 1
  subplot(1,3,3);
  imagesc(mod(CL21.chrn,2));
  ax(3)=gca;
  set(gca,'XTick',[]);
  set(gca,'YTick',chrnpos,'YTickLabel',CL21.chr(chrnpos));
end

ax=ax(~isnan(ax));

if length(ax)>1
  linkprop(ax,'CLim');
  caxis([-1 1]);
  linkaxes(ax);
end




