function run_somcall(sample,params)

if ~exist('params','var'), params=[]; end
params=impose_default_value(params,'lane_blacklist_file',[]);
params=impose_default_value(params,'lanetable_file',[]);

t_boomdir = ['/xchip/tcga_scratch/ng/' sample '/wgs/boom/tumor.boom'];
n_boomdir = ['/xchip/tcga_scratch/ng/' sample '/wgs/boom/normal.boom'];
outdir = ['/xchip/tcga_scratch/ng/' sample '/wgs/somcall'];

scatter_somcall(t_boomdir,n_boomdir,outdir,sample,params);

