function [lnk,idx] = cl2_linkage(clustering_params,similarity,D)

maxD=max(D(:));
if similarity
  D = maxD - D;
end

if ischar(clustering_params)
  cl_type=clustering_params;
else
  cl_type=clustering_params.method;
end

switch lower(cl_type)
 case {'single','average','complete','centroid','ward','weighted'}
  if (0) % exist('flinkage')==3
    D2=D;
    [lnk,idx] = flinkage(D2, cl_type );
  else
    [lnk,idx] = z2lnk( linkage( mat2pdist( D ), cl_type ) );
  end
 case 'consensus'
  [lnk,idx] = consensus_clustering(clustering_params,similarity,D);
 case 'spc'  
 case 'mfspc'
 case 'sspc'
  [Jij,a]=create_spc_Jij(dat,clustering_params.K);
 otherwise
  disp('Unknown method');
end

if similarity
  lnk(:,5)=maxD - lnk(:,5);
end




