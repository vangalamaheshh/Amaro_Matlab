function lh=cum_dist(x)

for i=1:length(x)
  x{i}=as_column(x{i});
end

all=cat(1,x{:});
rmin=min(all);
rmax=max(all);

clf;
cm=colormap;
lh=[];
for i=1:length(x)
  sx=sort(x{i});
  [u,ui,uj]=unique(sx);
  u=[ rmin; u; rmax];
  ui=[ 0; ui; length(sx)]/length(sx);
  lnx=[];
  lny=[];
  for j=1:(length(u)-1);
    lnx=[lnx u(j) u(j+1) u(j+1)];
    lny=[lny ui(j) ui(j) ui(j+1)];
  end
  lh(i)=line(lnx,lny,'Color',cm(i,:));
end
axis([rmin rmax 0 1]);


