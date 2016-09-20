function lhs=show_dend_graph(g,val,mbr,c,N,nan_val,lw,hf)
% g - graph
% val - value per node
% mbr - members of each node
% c - color of each node
% N - number of points
% nan_val - what to do with NaN
% p - patches coordinates

if ~exist('hf','var')
  hf=-ones(size(c));
end

val(find(isnan(val)))=nan_val;
low=min(val);
high=max(val);
np=find(sum(spones(g),2)>0);

nonempty_rows=find(sum(spones(g),2));
dum=nonempty_rows(1);
if val(find(g(dum,:))) > val(dum)
  flip_patch=0;
else
  flip_patch=1;
end
set(gca,'UserData','Dendrogram');

global dend_patch_width
if ~isempty(dend_patch_width)
   xsz=N*dend_patch_width;
else
   xsz=N/64;
end
hgt=(high-low)/64;
axis([0.5 N+0.5 low-hgt high+hgt]);
if (isempty(c))
  c = ones(size(mbr,1),1);
end;

if ( min(c) == max(c) )
  caxis([min(c)-1 min(c)+1])
else
  caxis([min(c) max(c)]);
end
nn=size(g,1);
[x,y,e]=find(g);
lhs=[];
for i=1:length(x)
    x1=mean(mbr(x(i),:));
    x2=mean(mbr(y(i),:));
    y1=val(x(i));
    y2=val(y(i));
    cury=min( hgt, abs(y2-y1)*0.8);
    curyf=cury;
    if ( hf(x(i)) >= 0 )
      curyf=curyf*hf(x(i));
    end
    curx=min( 0.8*(mbr(x(i),2)-mbr(x(i),1)+1), xsz )/2;
%    curx=(mbr(x(i),2)-mbr(x(i),1)+1)/2;
    if flip_patch
      cury=-cury;
      curyf=-curyf;
    end
    px1=x1-curx; px2=x1+curx;
    py1=y1+cury-curyf; py2=y1+cury;
    l1=line( [ x1 x1 x2 ],[y1 y2 y2]);
    if exist('lw','var') && ~isempty(lw) 
      set(l1,'LineWidth',lw);
    end
    set(l1,'Zdata',[0 0 0]); % was [ 0 0 ]
    set(l1,'EraseMode','background');
%    l2=line( [ x1 x2 ],[y2 y2 ]); 
%    if exist('lw','var') && ~isempty(lw)
%      set(l2,'LineWidth',lw);
%    end
%    set(l2,'Zdata',[0 0]);
%    set(l2,'EraseMode','background');
    lhs=[lhs l1];
    if c(x(i)) >= 0 
      h=patch( [px1 px2 px2 px1 px1], ...
	       [py1 py1 py2 py2 py1], c(x(i)));
      set(h,'UserData',[ x(i) val(x(i)) mbr(x(i),:) c(x(i)) ]);
      set(h,'Zdata',[ 1 1 1 1 1]);
    end
end
% zoom on
