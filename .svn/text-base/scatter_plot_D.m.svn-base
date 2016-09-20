function scatter_plot_D(X,D,color_sid,marker_sid)

hold on
if exist('color_sid','var') && ~isempty(color_sid)
  if ~exist('marker_sid','var') || isempty(marker_sid)
    marker_sid=color_sid;
  end
  [D1,markcol]=crossprod_types(D,marker_sid,color_sid);
  ut=unique(D1.supdat(markcol,:));

  for i=1:length(ut)
    idx=find(D1.supdat(markcol,:)==ut(i));
    if size(X,2)==2
      ph{i}=plot(X(idx,1),X(idx,2),'.');
    else
      ph{i}=plot3(X(idx,1),X(idx,2),X(idx,3),'.');
    end  
    set(ph{i},'MarkerSize',D.supmark(marker_sid).marker_size(D.supdat(marker_sid,idx(1))),...
              'LineWidth',D.supmark(marker_sid).marker_size(D.supdat(marker_sid,idx(1)))*0.1,...
              'Marker',D.supmark(marker_sid).marker{D.supdat(marker_sid,idx(1))},...
              'MarkerFaceColor',D.supmark(color_sid).colormap(D.supdat(color_sid,idx(1)),:),...
              'MarkerEdgeColor',D.supmark(color_sid).colormap(D.supdat(color_sid,idx(1)),:));
  end
else
  if size(X,2)==2
    ph=plot(X(:,1),X(:,2),'.');
  else
    ph=plot3(X(:,1),X(:,2),X(:,3),'.');
  end  
end
