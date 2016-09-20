function [Dord,dend]=one_way_clustering(D,rowcol,params)
%  [Dord,dend]=one_way_clustering(D,rowcol,params)
%  D - A data object
%  rowcol (string) - input 'col' or 'row' to indicate which dimension to cluster
%  params (struct) - params.cluster: clustering method (string) e.g. 'average', 'complete' ...
%                    params.dist: distance measure as needed for function dist.m
%                    params.distmat (optional): use this matrix as a distace matrix
%                    params.is_sim  (optional): is distmat a similarity measure? (highest is closest)
%                    params.preproc (optional): run preprocess_D with this parameter before clustering
%
%  Dord - ordered version of D
%  dend - dend.idx: reordering vector for D
%         dend.lnk: (n-1)x5 (n=no. of objected) [ L_from:L_to R_from:R_to val ]

if ~exist('params','var')
  params=struct('cluster','average','dist','correlation');
end

D1=D;
verbose('Preprocessing ....');
if isfield(params,'preproc')
  D1=preprocess_D(D,params.preproc);
end

is_sim=0;
if isfield(params,'distmat')
  dist=params.distmat;
  is_sim=params.is_sim;
else
  switch rowcol(1:3)
   case {'col','sam','con','exp'}
    verbose(['Calculating distance based on ' num2str(size(D1.dat,1)) ' features']);
    dist=dna_distance(D1.dat',params.dist);
   case {'row','gen','mir'}
    verbose(['Calculating distance based on ' num2str(size(D1.dat,2)) ' features']);
    dist=dna_distance(D1.dat,params.dist);
   otherwise
    error('no match');
  end
end
verbose(['About to cluster ' num2str(size(dist,1)) ' objects']);
[dend.lnk,dend.idx]=cl2_linkage(params.cluster,is_sim,dist);
if isfield(params,'order')
  [dend.lnk,dend.idx]=order_dend(params.order,dend.lnk,dend.idx,dist);
end
[dend.g,dend.v,dend.m]=lnk_to_graph(dend.lnk,length(dend.idx));
[dum,dend.revidx]=sort(dend.idx);
dend.z=lnk2z(dend.lnk,dend.idx);
dend.params=params;
Dord=reorder_D(D,rowcol,dend.idx);
