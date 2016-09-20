function [lnk,idx,Cmean,C,run]=consensus_clustering(p,similarity,D)

% repeat n_iter times w/ or w/o lsf
% introduce noise: subsample objects, (subsample features - can't do
% it - no access to features), reruns with different random seed
% count coclustered at level x according to nclusters or height at dendrogram
% cluster the outcome

N=size(D,1);
for i=1:p.n_iter
  switch p.noise.method
   case 'subsample_objects'
    nobjs=round(p.noise.frac*N);
    if nobjs<2
      error(['Need at least 2 objects']);
    end
    r=randperm(N);
    cobjs=r(1:nobjs);
    run(i).objs=cobjs;
    [run(i).lnk,run(i).idx]=cl2_linkage(p.repeated_clustering_params, ...
                                        similarity,D(cobjs,cobjs));
   case 'random_seed'
  end
end

% count co-clustered
C=NaN*ones(N,N,p.n_iter);
for i=1:p.n_iter
  switch p.combine.method
   case 'nclusters'
    coph=lnk2cophenet(run(i).lnk,1);
    p.combine.is_similarity=1;
   case 'height'
    coph=lnk2cophenet(run(i).lnk,0);    
    p.combine.is_similarity=similarity;
   case 'k'
    coph=lnk2cophenet(run(i).lnk,1);
    coph=(coph>=p.combine.k);
    p.combine.is_similarity=1;    
  end
  [tmp,revidx]=sort(run(i).idx);
  coph=coph(revidx,revidx);
  run(i).coph=coph;
  C(run(i).objs,run(i).objs,i)=coph;
end

Cmean=nanmean(C,3);
[lnk,idx]=cl2_linkage(p.post_clustering_params,p.combine.is_similarity,Cmean);

