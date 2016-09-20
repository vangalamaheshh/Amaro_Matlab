# filter miGCM_218.gct
gp_cluster miGCM_218.filt.gct miGCM_218.filt.sord.gct miGCM_218.filt.sdend.odf 1 1 2 1    

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

gp_cluster miGCM_218.filt.sord.gct miGCM_218.filt.sord.gord.gct miGCM_218.filt.gdend.odf 1 1 2 0

gp_clustering_figure -o test.png -d miGCM_218.filt.sord.gord.gct -dp 1 ...
    -p Tissue Type -c miGCM_218.tt.cls ...
    -cc color_table.txt -gd miGCM_218.filt.gdend.odf -gdr 1 -sd ...
    miGCM_218.filt.sdend.odf -sdr 0 ...
    -sf 5 -gf 10 -sx 0.3 0.2 1 -sy 0.3 0.1 0.1 0.4 0.05 -l 1 


-dp data_preprocessing=enum_param(a.dp{1},...
                                  {0,'none'; 1, ...
                        'row_centerd_and_normalize';2,'col_centered_and_normalize'});
-p phenotype 
-gdr reduce_denrogram    gdend_reduced=enum_param(a.gdr{1},{0,'all';1,'reduced_big'});


-sf sample font size
-gf gene font size
-sx size x
-sy size y
-l landscape


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_features TN.new.gct TN.filt.gct 7.25
gp_cluster TN.filt.gct TN.filt.sord.euc.gct TN.filt.sdend.euc.odf 1 1 0 1    
gp_cluster TN.filt.sord.euc.gct TN.filt.sord.euc.gord.gct TN.filt.gdend.odf 1 1 2 0

gp_clustering_figure -o test2.png -d TN.filt.sord.euc.gord.gct -dp 1 ...
    -p -c TN.new.cls ...
    -cc color_table.txt -gd TN.filt.gdend.odf -gdr 1 -sd ...
    TN.filt.sdend.euc.odf -sdr 0 ...
    -sf 5 -gf 10 -sx 0.3 0.2 1 -sy 0.3 0.1 0.1 0.4 0.05 -l 1 
 
