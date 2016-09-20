function display_confmat(D,supid,confmat,should_cluster,idx)
[typeacc,typedesc,Y,range,non_empty]=decollapse_supdat(D,supid,1:max(D.supdat(supid,:)),0);
X=make_D(confmat,Y.supacc(range,:),[],Y.supacc(range,:));
X.dat=X.dat./repmat(sum(X.dat,2),1,size(X.dat,2));
X.sdesc=strvcat(X.sdesc);
X.gacc=strvcat(X.gacc);
if exist('should_cluster','var') && should_cluster
  [Xord,dend]=one_way_clustering(X,'col',struct('cluster','average','distmat',1-X.dat,'is_sim',0));
  Xord=reorder_D_rows(Xord,Xord.origidx);
else
  if exist('idx','var') && ~isempty(idx)
    Xord=reorder_D_rows(reorder_D_cols(X,idx),idx);
  else
    Xord=X;
  end
  dend=[];
end
disp_p=struct('x',struct('sizes',[0.35 0.35 1],'gaps',[1 0.1 0.1 1],'border',0.1),...
              'y',struct('sizes',[0.35 1],'gaps',[1 0.1 1],'border',0.1),...
              'items',{{...
                  {2,1,'gdend','all',0.5},...
                  {1,3,'sdesc','horizontal',10,[],0,1,{'HorizontalAlignment','Left'}},...
                  {2,2,'gacc','vertical',10},...
                  {2,3,'consensusmat',[0 max(Xord.dat(:))]}}});
display_D(Xord,dend,[],disp_p);
