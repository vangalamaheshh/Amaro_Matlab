function plot_sa(tn,num)

fname = ['/xchip/tcga_scratch/lawrence/mm/0309/wgs/' tn '_dR2/sa1/results_batch' ...
  num2str(num) '.txt'];

T = load_matrix(fname);
scatter(T(:,7)./T(:,3),T(:,8));
