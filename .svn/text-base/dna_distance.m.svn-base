function d=dna_distance(vecs,dist_params,varargin)

if ischar(dist_params)
  dist_type=dist_params;
else
  dist_type=dist_params.method;
end

try
  switch dist_type
   case 'cn_euclid'
     d=dna_distance(dna_norm(vecs')','euclid');
   case 'cn_correlation'
     d=dna_distance(dna_norm(vecs')','correlation');    
   case 'nanjaccard'
    d=1-jaccard(vecs);
   otherwise
    d=squareform(pdist(vecs,dist_type,varargin{:}));
  end
catch
  disp('dna_distance: out of memory ... ');
  if strcmp(dist_type,'euclid')
    disp('...using dna_dist2');
    d=dna_dist2(vecs);
  else
    d=zeros(size(vecs,1),size(vecs,1));
    for i=1:size(vecs,1)
      d(i,:)=dist(vecs(i,:),vecs,dist_type);
    end
  end
end
