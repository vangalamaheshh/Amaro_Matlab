function chromosome_figures(fname,C,cyto,cluster_params,disp_params)

if ~isfield(C,'cyton')
  C=add_cyto(C,cyto);
end

if ~exist('cluster_params','var') || isempty(cluster_params)
  cluster_params=struct('cluster','average','dist','cn_euclid');
end

if ~exist('disp_params','var') || isempty(disp_params)
%   disp_params={'snpcyto',struct('y',struct('sizes',[ 0.5 0.025],'gaps',[1 0.2 1],'border',0.1),...
%                             'items',{{{1,4,'dataorig',[-1 1]},...
%                       {2,4,'colorbar','horizontal','for',1,4,4,6}}})};
  disp_params={'snpcytosup',struct('y',struct('sizes',[ 0.025 0.5 0.025],'gaps',[1 0.2 1],'border',0.1),...
                            'items',{{{2,4,'dataorig',[-1 1]},...
                      {3,4,'colorbar','horizontal','for',1,4,4,6}}})};
end

close all
for i=1:max(C.chrn)
  figure(i);
  in_chr=find(C.chrn==i);
  if ~isempty(in_chr)
    X=reorder_D_rows(C,in_chr);
    svdX=svd_D(X,'cols',1);
    [tmp,idx]=sort(svdX.dat(1,:));
    Xord=reorder_D_cols(X,idx);
%    [Xord,Xord.sdend]=one_way_clustering(X,'cols',cluster_params);
    display_D(Xord,[],[],disp_params);
    if ~isempty(fname)
      print_D([ fname '.chr' num2chromosome(i)],{{'fig'},{'pdf'},{'png','-r180'}});
      close(i);
    end
  end
end
