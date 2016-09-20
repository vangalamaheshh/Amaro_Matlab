function [ph,lh,D,ln,sh]=plot_D(D,rc,ids,color_supid,shape_supid,add_lines_to_plane,add_shadow)

ndim=length(ids);
if length(ids)<2 || length(ids)>3
  error('number of ids not supported');
end

if is_col(rc)
  dat=D.dat';
else
  dat=D.dat;
end

no_supid=0;
if ~exist('color_supid','var') || isempty(color_supid)
  s=ones(size(dat,1),1);
  u=1;
  marker='.';
  markersize=10;
  markercolor=[ 0 0 1];
  no_supid=1;
elseif ~exist('shape_supid','var') || isempty(shape_supid)
  s=D.supdat(color_supid,:);
  u=unique(s);
  u=u(find(~isnan(u) & u~=0));
  shape_supid=color_supid;
else
  [D,colorshape_supid]=crossprod_types(D,color_supid,shape_supid);
  s=D.supdat(colorshape_supid,:);
  u=unique(s);
  u=u(find(~isnan(u) & u~=0));
end

ln={};
sh={};
if ndim==2
  for i=1:length(u)
    idx=find(s==u(i));
    ph{i}=plot(dat(idx,ids(1)),dat(idx,ids(2)),'.'); hold on
    if no_supid
      set(ph{i},'MarkerSize',markersize,...
                'Marker',marker,...
                'MarkerFaceColor',markercolor,...
                'MarkerEdgeColor',markercolor);
    else
      set(ph{i},'MarkerSize',D.supmark(shape_supid).marker_size(D.supdat(shape_supid,idx(1))),...
                'Marker',D.supmark(shape_supid).marker{D.supdat(shape_supid,idx(1))},...
                'MarkerFaceColor',D.supmark(color_supid).colormap(D.supdat(color_supid,idx(1)),:),...
                'MarkerEdgeColor',D.supmark(color_supid).colormap(D.supdat(color_supid,idx(1)),:));
    end
  end
  if is_col(rc)
    if ischar(D.gacc)
      D.gacc=cellstr(D.gacc);
    end
    lh(1)=xlabel(D.gacc{ids(1)});
    lh(2)=ylabel(D.gacc{ids(2)});
  else
    if ischar(D.sdesc)
      D.sdesc=cellstr(D.sdesc);
    end
    lh(1)=xlabel(D.sdesc{ids(1)});
    lh(2)=ylabel(D.sdesc{ids(2)});
  end
  set(lh,'Interpreter','None');
  axis equal
  axis fill
elseif ndim==3
  a=[min(dat(find(ismember(s,u)),ids),[],1); max(dat(find(ismember(s,u)),ids),[],1)];
  axis(a(:)');
  enlarge_axis(0.05);
  box on
  axis equal
  axis vis3d
  ax=axis;
  if exist('add_lines_to_plane','var')
    for i=1:length(u)
      idx=find(s==u(i));
      for j=1:length(idx)
        tmp=[dat(idx(j),ids) ;dat(idx(j),ids)];
        tmp(1,add_lines_to_plane.axis)=ax(add_lines_to_plane.axis*2-1);
        ln{i}(j)=line(tmp(:,1),tmp(:,2),tmp(:,3)); hold on;
        set(ln{i}(j),'Color',add_lines_to_plane.color,'EraseMode','background');
      end
    end
  end
  for i=1:length(u)
    idx=find(s==u(i));
    ph{i}=plot3(dat(idx,ids(1)),dat(idx,ids(2)),dat(idx,ids(3)),'.'); hold on
    if no_supid
      set(ph{i},'MarkerSize',markersize,...
                'Marker',marker,...
                'MarkerFaceColor',markercolor,...
                'MarkerEdgeColor',markercolor);
    else
      if D.supmark(shape_supid).marker{D.supdat(shape_supid,idx(1))}=='S'
        make_spheres(ph{i},D.supmark(shape_supid).marker_size(D.supdat(shape_supid,idx(1))));
      else
        set(ph{i},'MarkerSize',D.supmark(shape_supid).marker_size(D.supdat(shape_supid,idx(1))),...
                  'Marker',D.supmark(shape_supid).marker{D.supdat(shape_supid,idx(1))},...
                  'MarkerFaceColor',D.supmark(color_supid).colormap(D.supdat(color_supid,idx(1)),:),...
                  'MarkerEdgeColor',D.supmark(color_supid).colormap(D.supdat(color_supid,idx(1)),:));
      end
    end
    if exist('add_shadow','var')  && exist('add_lines_to_plane','var')
      tmp=dat(idx,ids);
      tmp(:,add_lines_to_plane.axis)=ax(add_lines_to_plane.axis*2-1);
      sh{i}=plot3(tmp(:,1),tmp(:,2),tmp(:,3),'.'); hold on
      if no_supid
        set(sh{i},'MarkerSize',add_shadow.markersize,...
                  'Marker',add_shadow.marker,...
                  'MarkerFaceColor',markercolor,...
                  'MarkerEdgeColor',markercolor);
      else
        set(sh{i},'MarkerSize',add_shadow.markersize,...
                  'Marker',add_shadow.marker,...
                  'MarkerFaceColor',D.supmark(color_supid).colormap(D.supdat(color_supid,idx(1)),:),...
                  'MarkerEdgeColor',D.supmark(color_supid).colormap(D.supdat(color_supid,idx(1)),:));      
      end
    end
  end
  if is_col(rc)
    if ischar(D.gacc)
      D.gacc=cellstr(D.gacc);
    end
    lh(1)=xlabel(D.gacc{ids(1)});
    lh(2)=ylabel(D.gacc{ids(2)});
    lh(3)=zlabel(D.gacc{ids(3)});
  else
    if ischar(D.sdesc)
      D.sdesc=cellstr(D.sdesc);
    end
    lh(1)=xlabel(D.sdesc{ids(1)});
    lh(2)=ylabel(D.sdesc{ids(2)});
    lh(3)=zlabel(D.sdesc{ids(3)});
  end
  set(lh,'Interpreter','None');
end
hold off

