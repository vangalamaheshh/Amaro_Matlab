function display_D_dist(D,rc,params)

if isfield(params,'preproc')
  D1=preprocess_D(D,params.preproc);
else
  D1=D;
end

if is_col(rc)
  dist=dna_distance(D1.dat',params.dist);
else
  dist=dna_distance(D1.dat,params.dist);
end

imagesc(dist);
% colormap default
