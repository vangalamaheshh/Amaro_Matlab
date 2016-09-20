function mdsD=mds_D(D,rc,k,dist_type,mds_type)

if ~exist('dist_type','var') || isempty(dist_type)
  dist_type='euclidean';
end

if ~exist('mds_type','var')|| isempty(mds_type)
  mds_type.method='metricstress';
  mds_type.start='random';
end

if is_col(rc)
  distmat=dist(D.dat',[],dist_type);
%   keyboard
  [m,stress,disparities]=mdscale(distmat,k,'Criterion',mds_type.method,'start',mds_type.start,'Options',mds_type.opts);
  mdsD=reorder_D_rows(D,[]);
  mdsD.dat=m';
  mdsD.gacc=cellstr([repmat('MDS #',k,1) num2str((1:k)','%d')]);
  mdsD.gdesc=mdsD.gacc;
  mdsD.mds.m=m;
  mdsD.mds.stress=stress;
  mdsD.mds.disparities=disparities;
else
  [m,stress,disparities]=mdscale(dist(D.dat,[],dist_type),k,'Criterion',mds_type.method,'start',mds_type.start,'Options',mds_type.opts);
  mdsD=reorder_D_cols(D,[]);
  mdsD.dat=prj;
  mdsD.sdesc=cellstr([repmat('MDS #',k,1) num2str((1:k)','%d')]);
  mdsD.mds.m=m;
  mdsD.mds.stress=stress;
  mdsD.mds.disparities=disparities;
end
