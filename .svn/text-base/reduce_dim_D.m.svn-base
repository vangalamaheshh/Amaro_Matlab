function R=reduce_dim_D(D,rc,k,red_dim_type)

if ischar(red_dim_type)
  tmp.method=red_dim_type;
  red_dim_type=tmp;
end

switch red_dim_type.method
 case 'PCA'
  R=pca_D(D,rc,k);
 case 'SVD'
  R=svd_D(D,rc,k);
 case 'MDS'
  R=mds_D(D,rc,k);
 otherwise
  error('no such type');
end

