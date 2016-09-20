function y=ctr_norm(x,ctr_type,norm_type)

if ischar(ctr_type)
  tmp.method=ctr_type;
  ctr_type=tmp;
end

if ischar(norm_type)
  tmp.method=norm_type;
  norm_type=tmp;
end

n=size(x,2);

switch ctr_type.method
 case 'mean'
  x1=x-repmat(mean(x,2),1,n);
 case 'median'
  x1=x-repmat(median(x,2),1,n);
 case 'prctile'
  x1=x-repmat(prctile(x,ctr_type.prctile,2),1,n);
 case 'min'
  x1=x-repmat(min(x,[],2),1,n);
 otherwise
  error('no such type');
end

switch norm_type.method
 case 'std'
  y=x1./repmat(std(x1,0,2),1,n);
 case 'unitvec'
  y=x1./repmat(std(x1,0,2)*sqrt(n-1),1,n);
 case 'mad'
  y=x1./repmat(mad(x1,1,2),1,n);
 case 'range'
  y=x1./repmat(range(x1,2),1,n);
 otherwise
  error('no such type');
end
