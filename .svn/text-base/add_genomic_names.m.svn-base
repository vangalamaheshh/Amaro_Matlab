function h=add_genomic_names(D,cyto,axis,fs)
N=40;
h=[];
if lower(axis)=='x'
  set(gca,'XTick',[]);
  ax=axis;
  for i=1:round(length(D.cyto)/N):length(D.cyto)
    text(i,0,cyto(D.cyto(i)).name,'FontSize',fs,'Rotation',270,'HorizontalAlignment','left');
  end
else
  set(gca,'YTick',[]);
  ax=axis;
  for i=1:round(length(D.cyto)/N):length(D.cyto)
    hi=text(0,i,cyto(D.cyto(i)).name,'FontSize',fs, ...
            'HorizontalAlignment','right');
    h=[h; hi];
  end
end

