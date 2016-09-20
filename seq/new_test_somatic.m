% simple test

n_boomdir = '/xchip/tcga_scratch/lawrence/gbm/0188/wgs/normal.boom';
t_boomdir = '/xchip/tcga_scratch/lawrence/gbm/0188/wgs/tumor.boom';
chr = 17;
start = 7e6;
stop = 8e6;
params = [];
M = somcall(t_boomdir,n_boomdir,chr,start,stop,params);
T = annotate_M_calls(M);


% call whole exome

sample = 'ov/1319/wgs';
n_boomdir = ['/xchip/tcga_scratch/lawrence/' sample '/normal.boom'];
t_boomdir = ['/xchip/tcga_scratch/lawrence/' sample '/tumor.boom'];
outdir = ['/xchip/tcga_scratch/lawrence/' sample '/somcall'];
params = [];
%params.method = 'pv';
scatter_somcall(t_boomdir,n_boomdir,outdir,params);












% debugging



%%% debug failed chr
n_boomdir = ['/xchip/tcga_scratch/lawrence/' sample '/normal.boom'];
t_boomdir = ['/xchip/tcga_scratch/lawrence/' sample '/tumor.boom'];
outdir = ['/xchip/tcga_scratch/lawrence/' sample '/somcall'];
params = [];
%params.method = 'pv';
chr = 15;
scatter_somcall_chr(t_boomdir,n_boomdir,outdir,chr,params);



%%% debug failed chunk

% load exome
b=read_table(['/xchip/tcga_scratch/gadgetz/mutation_analysis/' ...
              'whole_exome_refseq_coding.targets.interval_list.no_header'],...
             '%s%f%f%*s%*s',char(9),0);
X.chrn=convert_chr(b.dat{1});
X.start=b.dat{2};
X.end=b.dat{3};
chrlen=load_chrlen;

n_boomdir = '/xchip/tcga_scratch/lawrence/gbm/0188/wgs/normal.boom';
t_boomdir = '/xchip/tcga_scratch/lawrence/gbm/0188/wgs/tumor.boom';
outdir = '/xchip/tcga_scratch/lawrence/gbm/0188/wgs/somcall';
chr = 5;
params = [];
%params.method = 'pv';
X2 = reorder_struct(X,X.chrn==chr & X.end<chrlen(chr));
first = 1;
last=min(slength(X2),first+999);
fprintf('\n\nTARGETS %d-%d of %d\n\n',first,last,slength(X2));
M = somcall(t_boomdir,n_boomdir,X2.chrn(first:last),X2.start(first:last),X2.end(first:last),params);


