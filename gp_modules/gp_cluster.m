function gp_cluster(gctfile,out_gctfile,out_dendfile,...
                    preproc_type,clustering_type,distance_type,rc,use_treeview_format)

params.preproc=enum_param(preproc_type,...
                          {0 'none'; ...
                    1 'row_center_and_normalize'; ...
                    2 'col_center_and_normalize'});

params.cluster=enum_param(clustering_type,...
                          {0,'single'; 1,'average';2,'complete'; ...
                    3,'weighted';4,'centroid';5,'median';6,'ward'});

params.dist=enum_param(distance_type,...
                       {0,'euclidean';1,'cosine';2,'correlation'});

rowcol=enum_param(rc,{0,'row';1,'col'});

if ~exist('use_treeview_format')
  use_treeview_format=0;
end

% read gctfile
D=read_mit_gct_file(gctfile);

[Dord,dend]=one_way_clustering(D,rowcol,params);

if use_treeview_format
  if dend.lnk(end,5)>dend.lnk(1,5)
    dend.lnk(:,5)=2-dend.lnk(:,5);
  end
  write_mit_gct_file(out_gctfile,Dord);
  if is_col(rowcol) %col
    [lnk,idx]=blank_lnk_idx(size(D.dat,1));
    write_cdt_gtr_atr(out_dendfile,out_gctfile,lnk,idx,dend.lnk,dend.idx);
  else %row
    [lnk,idx]=blank_lnk_idx(size(D.dat,2));
    write_cdt_gtr_atr(out_dendfile,out_gctfile,dend.lnk,dend.idx,lnk,idx);
  end    
else
  write_mit_gct_file(out_gctfile,Dord);
  write_mit_dend_file(out_dendfile,dend);
end
