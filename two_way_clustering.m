function [D2,gdend,sdend]=two_way_clustering(D,gparams,sparams,preproc_params)

if exist('preproc_params','var') && ~isempty(preproc_params)
  if ischar(preproc_params)
    tmp.method=preproc_params;
    preproc_params=tmp;
  end
  if strcmp(preproc_params.method,'rowcolnorm')
    D.dat=row_and_column_norm(D.dat);
  end
%    dist_type='euclid';
end

if ~exist('gparams','var') || isempty(gparams)
  gparams=struct('cluster','average','dist','correlation');
end
if ~exist('sparams','var') || isempty(sparams)
  sparams=struct('cluster','average','dist','correlation');
end

[D1,sdend]=one_way_clustering(D,'samples',sparams);
[D2,gdend]=one_way_clustering(D1,'genes',gparams);

if nargout==1
  D2.sdend=sdend;
  D2.gdend=gdend;
end





