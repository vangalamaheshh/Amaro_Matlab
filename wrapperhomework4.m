% Summarizes k-means and kpp output
cd /Users/sbiswas/GitHub/cs181-s16-homework1-surgebiswas/homework4/kmeans_output
!rm agg*
aggregate_plots('kmeans_mean', true, false, 1);
 
for i = 0 : 9
    aggregate_plots(['kmeans_rep_cluster_', num2str(i), '_'], true, false, 1);
end
 
aggregate_plots('agg_kmeans_rep_cluster', false, true, 0.1);